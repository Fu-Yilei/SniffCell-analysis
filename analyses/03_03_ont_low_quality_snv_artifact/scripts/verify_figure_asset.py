#!/usr/bin/env python3
"""Verify the exact manuscript Supplementary Figure 12 image."""

from hashlib import sha256
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ASSET = HERE / "results" / "supplementary_figure_12.png"
EXPECTED = "7464ab67950ebd5943f1c2c083716412744f9b40d45cb3cc6e4c600fe7c71aa7"

digest = sha256(ASSET.read_bytes()).hexdigest()
if digest != EXPECTED:
    raise SystemExit(f"Checksum mismatch: expected {EXPECTED}, observed {digest}")
print(f"OK  {digest}  {ASSET.name}")
