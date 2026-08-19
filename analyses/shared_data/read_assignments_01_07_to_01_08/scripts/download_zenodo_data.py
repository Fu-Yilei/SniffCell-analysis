#!/usr/bin/env python3
"""Download, verify, and extract the shared SniffCell read assignments."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path, PurePosixPath

import requests


RECORD_ID = 22003085
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            result.update(chunk)
    return result.hexdigest()


def parse_args() -> argparse.Namespace:
    shared_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=shared_dir / "data")
    parser.add_argument("--keep-archive", action="store_true")
    return parser.parse_args()


def remote_files() -> dict[str, str]:
    response = requests.get(API_URL, timeout=60)
    response.raise_for_status()
    entries = response.json().get("files", {}).get("entries", {})
    return {
        name: entry.get("links", {}).get("content", "")
        for name, entry in entries.items()
    }


def validate_members(archive: Path, expected_root: str) -> None:
    listing = subprocess.run(
        ["tar", "-I", "zstd", "-tf", str(archive)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    if not listing:
        raise RuntimeError("Archive is empty")
    for member in listing:
        path = PurePosixPath(member)
        if path.is_absolute() or ".." in path.parts or path.parts[0] != expected_root:
            raise RuntimeError(f"Unsafe archive member: {member}")
    table_count = sum(name.endswith("/deconv_read_summary.tsv") for name in listing)
    if table_count != 22:
        raise RuntimeError(f"Expected 22 read-assignment tables, found {table_count}")


def main() -> int:
    args = parse_args()
    shared_dir = Path(__file__).resolve().parents[1]
    with (shared_dir / "archive_manifest.tsv").open(newline="") as handle:
        archive_row = next(csv.DictReader(handle, delimiter="\t"))

    args.data_dir.mkdir(parents=True, exist_ok=True)
    archive = args.data_dir / archive_row["zenodo_filename"]
    archive_ok = (
        archive.is_file()
        and archive.stat().st_size == int(archive_row["bytes"])
        and digest(archive) == archive_row["sha256"]
    )
    if not archive_ok:
        url = remote_files().get(archive.name)
        if not url:
            raise FileNotFoundError(
                f"Zenodo record is missing {archive.name!r}; the draft may not be published."
            )
        temporary = archive.with_suffix(archive.suffix + ".part")
        with requests.get(url, stream=True, timeout=(60, 7200)) as response:
            response.raise_for_status()
            with temporary.open("wb") as output:
                for chunk in response.iter_content(16 * 1024 * 1024):
                    if chunk:
                        output.write(chunk)
        if temporary.stat().st_size != int(archive_row["bytes"]) or digest(temporary) != archive_row["sha256"]:
            temporary.unlink(missing_ok=True)
            raise RuntimeError("Archive checksum verification failed")
        temporary.replace(archive)

    validate_members(archive, archive_row["root_folder"])
    subprocess.run(
        ["tar", "-I", "zstd", "-xf", str(archive), "-C", str(args.data_dir)],
        check=True,
    )
    print(f"verified and extracted 22 tables under {archive_row['root_folder']}/")
    if not args.keep_archive:
        archive.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
