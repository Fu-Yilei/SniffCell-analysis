#!/usr/bin/env python3
"""Download, extract, and verify the shared modkit-derived CpG pileups."""

from __future__ import annotations

import argparse
import csv
import hashlib
import tarfile
from pathlib import Path

import requests


RECORD_ID = 22003085
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"


def digest(path: Path, algorithm: str = "sha256") -> str:
    result = hashlib.new(algorithm)
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
    files = response.json().get("files", {})
    entries = files.get("entries", {}) if isinstance(files, dict) else {}
    return {
        name: entry.get("links", {}).get("content", "")
        for name, entry in entries.items()
    }


def safe_extract(archive: Path, destination: Path, expected_root: str) -> None:
    destination = destination.resolve()
    with tarfile.open(archive) as handle:
        members = handle.getmembers()
        roots = {Path(member.name).parts[0] for member in members if member.name}
        if roots != {expected_root}:
            raise RuntimeError(f"Unexpected archive roots: {sorted(roots)}")
        for member in members:
            target = (destination / member.name).resolve()
            if destination != target and destination not in target.parents:
                raise RuntimeError(f"Unsafe archive member: {member.name}")
        handle.extractall(destination)


def main() -> int:
    args = parse_args()
    shared_dir = Path(__file__).resolve().parents[1]
    with (shared_dir / "archive_manifest.tsv").open(newline="") as handle:
        archive_row = next(csv.DictReader(handle, delimiter="\t"))
    with (shared_dir / "file_manifest.tsv").open(newline="") as handle:
        file_rows = list(csv.DictReader(handle, delimiter="\t"))

    args.data_dir.mkdir(parents=True, exist_ok=True)
    archive = args.data_dir / archive_row["zenodo_filename"]
    if not (
        archive.is_file()
        and archive.stat().st_size == int(archive_row["bytes"])
        and digest(archive) == archive_row["sha256"]
    ):
        urls = remote_files()
        url = urls.get(archive.name)
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
    print(f"verified archive: {archive.name}")

    root = args.data_dir / archive_row["root_folder"]
    safe_extract(archive, args.data_dir, archive_row["root_folder"])
    for row in file_rows:
        path = root / row["relative_path"]
        if path.stat().st_size != int(row["bytes"]) or digest(path) != row["sha256"]:
            raise RuntimeError(f"File verification failed: {row['relative_path']}")
    print(f"verified {len(file_rows)} pileup files under {root.name}/")

    if not args.keep_archive:
        archive.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
