#!/usr/bin/env python3
"""Download section 1.2 files from Zenodo into the analysis folder layout."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd
import requests


RECORD_ID = 22003085
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=analysis_dir / "data/zenodo_manifest.tsv",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=analysis_dir / "data",
    )
    return parser.parse_args()


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def record_files(payload: dict) -> dict[str, str]:
    files = payload.get("files", {})
    entries = files.get("entries", {}) if isinstance(files, dict) else {}
    result: dict[str, str] = {}
    for name, entry in entries.items():
        links = entry.get("links", {})
        url = links.get("content") or links.get("self")
        if url:
            result[name] = url
    return result


def main() -> int:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    response = requests.get(API_URL, timeout=60)
    response.raise_for_status()
    remote = record_files(response.json())
    if not remote:
        raise RuntimeError(
            "No public files were returned. The Zenodo draft may not be published yet."
        )

    for row in manifest.itertuples(index=False):
        output = args.data_dir / row.relative_path
        expected_size = int(row.bytes)
        if (
            output.is_file()
            and output.stat().st_size == expected_size
            and sha256sum(output) == row.sha256
        ):
            print(f"verified {output.relative_to(args.data_dir)}")
            continue

        if row.zenodo_filename not in remote:
            raise FileNotFoundError(
                f"Zenodo record is missing {row.zenodo_filename!r}"
            )
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(output.suffix + ".part")
        with requests.get(remote[row.zenodo_filename], stream=True, timeout=(60, 7200)) as download:
            download.raise_for_status()
            digest = hashlib.sha256()
            size = 0
            with temporary.open("wb") as handle:
                for chunk in download.iter_content(chunk_size=16 * 1024 * 1024):
                    if not chunk:
                        continue
                    handle.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
        if size != expected_size or digest.hexdigest() != row.sha256:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(f"Checksum verification failed for {row.zenodo_filename}")
        temporary.replace(output)
        print(f"downloaded and verified {output.relative_to(args.data_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
