#!/usr/bin/env python3
"""Exploratory search for non-minimum-basis primitive sections of P48p.

First pass intentionally reports the exact archive schema and validates the
ambient Gram before selecting a search parameterization.
"""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "research/lower_bound_search_20261114/broad_dimensions_followup/ambient48.npz"
z = np.load(path, allow_pickle=False)
print("ARCHIVE", path, "bytes", path.stat().st_size)
for key in z.files:
    a = z[key]
    print("ARRAY", key, "shape", a.shape, "dtype", a.dtype,
          "min", np.min(a), "max", np.max(a))
