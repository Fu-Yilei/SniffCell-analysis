#!/usr/bin/env python3
"""Download the public per-locus supplement into a git-ignored private folder.

No source data are uploaded or committed. Use --expected-record-sha256 to reject
edits made since a reviewed snapshot. A hash mismatch is not a reason to alter
records: inspect the source revision and deliberately update the analysis.
"""
import argparse
from pathlib import Path
import tempfile
import urllib.request

from paper_source import read_candidates, records_sha256

URL = ("https://docs.google.com/spreadsheets/d/"
       "1o31UEyJkIofEj7v_beWsKueVgT5R-KwZPNxyhfV1ubw/export?format=csv&gid=1888934510")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("private_inputs/supplementary_table5.csv"))
    p.add_argument("--expected-record-sha256")
    a = p.parse_args()
    if "private_inputs" not in a.out.parts or a.out.suffix != ".csv":
        p.error("Output must be a CSV inside private_inputs/")
    if a.out.exists():
        p.error("Output exists; choose a new snapshot filename rather than overwriting it")
    with urllib.request.urlopen(URL, timeout=60) as response:
        content = response.read()
    with tempfile.TemporaryDirectory(prefix="table5-check-") as temporary:
        check = Path(temporary) / "source.csv"
        check.write_bytes(content)
        rows = read_candidates(check)
        digest = records_sha256(rows)
    if a.expected_record_sha256 and digest != a.expected_record_sha256:
        p.error("Table 5 changed: record checksum does not match the reviewed snapshot")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("xb") as handle:
        handle.write(content)
    print(f"Downloaded {len(rows)} candidate records; record SHA256 {digest}")


if __name__ == "__main__":
    main()
