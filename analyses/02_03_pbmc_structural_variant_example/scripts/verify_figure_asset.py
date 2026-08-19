#!/usr/bin/env python3
"""Verify the exact manually curated Supplementary Figure 14 asset."""

from hashlib import sha256
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ASSET = HERE / "results" / "supplementary_figure_14_source.png"
EXPECTED = "6dd6304d1681165559124a733a8ab2ee7ae0ad70e320905ebeedd0b8db3f7880"

digest = sha256(ASSET.read_bytes()).hexdigest()
if digest != EXPECTED:
    raise SystemExit(f"Checksum mismatch: expected {EXPECTED}, observed {digest}")
print(f"OK  {digest}  {ASSET.name}")
