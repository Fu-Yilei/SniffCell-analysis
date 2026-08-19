#!/usr/bin/env python3
"""Download and extract the Section 4 Zenodo archive."""

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path

RECORD_ID = "22003085"
ARCHIVE = "multi_region_brain_tr_figure_data.tar.zst"
SHA256 = "cad2a8d1eba03da99bd127b15f4e9e793d030173b1920e9d27e14715e29d6d99"
MEMBER = "multi_region_brain_tr_figure_data/pja2_read_lengths/pja2_read_lengths.tsv.gz"
ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--token-file", type=Path, help="Optional token for an unpublished draft")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data")
    args = parser.parse_args()
    headers = {}
    if args.token_file:
        headers["Authorization"] = f"Bearer {args.token_file.read_text().strip()}"
        api = f"https://zenodo.org/api/deposit/depositions/{RECORD_ID}"
    else:
        api = f"https://zenodo.org/api/records/{RECORD_ID}"
    with urllib.request.urlopen(urllib.request.Request(api, headers=headers), timeout=60) as fh:
        record = json.load(fh)
    entry = next(f for f in record["files"] if (f.get("key") or f.get("filename")) == ARCHIVE)
    url = entry["links"].get("content") or entry["links"].get("download")
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as fh:
        payload = fh.read()
    if hashlib.sha256(payload).hexdigest() != SHA256:
        raise ValueError("Downloaded archive SHA256 does not match the repository manifest.")
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / ARCHIVE
        archive.write_bytes(payload)
        zstd = subprocess.Popen(["zstd", "-dc", archive], stdout=subprocess.PIPE)
        subprocess.run(["tar", "-xf", "-", "-C", tmp], stdin=zstd.stdout, check=True)
        if zstd.stdout:
            zstd.stdout.close()
        if zstd.wait() != 0:
            raise RuntimeError("zstd extraction failed")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        destination = args.output_dir / Path(MEMBER).name
        shutil.copy2(Path(tmp) / MEMBER, destination)
    print(destination)


if __name__ == "__main__":
    main()
