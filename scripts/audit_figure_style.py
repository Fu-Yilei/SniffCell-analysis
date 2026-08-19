#!/usr/bin/env python3
"""Reject figure-level headers and footer notes in analysis plot scripts."""

from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PATTERNS = {
    "figure-level title": re.compile(r"\.(?:suptitle)\s*\(|\bplt\.title\s*\("),
    "figure-level footer text": re.compile(
        r"\b(?:fig|figure)\.text\s*\(|\bplt\.figtext\s*\("
    ),
}


def main() -> int:
    violations: list[tuple[Path, int, str, str]] = []
    for path in sorted((ROOT / "analyses").glob("*/scripts/*.py")):
        for line_number, line in enumerate(path.read_text().splitlines(), start=1):
            for label, pattern in PATTERNS.items():
                if pattern.search(line):
                    violations.append((path, line_number, label, line.strip()))

    if violations:
        for path, line_number, label, line in violations:
            print(f"FAIL {path.relative_to(ROOT)}:{line_number}: {label}: {line}")
        return 1

    print(
        "PASS: analysis scripts contain no figure-level top titles or bottom "
        "caption/footer text."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
