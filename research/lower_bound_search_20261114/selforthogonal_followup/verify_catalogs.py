"""Validate all downloaded generator matrices; completeness is a source claim."""

from pathlib import Path
import hashlib
import json
import numpy as np
from finite import ternary_code

P = Path(__file__).resolve().parent
data = np.load(P / "inputs/catalogs.npz")
rows = []
for n, t, expected in [
    (15, 7, 12),
    (16, 8, 7),
    (17, 8, 23),
    (18, 8, 160),
    (19, 9, 56),
    (20, 10, 24),
]:
    codes = data[str(n)]
    assert codes.shape == (expected, t, n)
    for G in codes:
        assert np.all((0 <= G) & (G < 3)) and not (G @ G.T % 3).any()
        C = ternary_code(G)
        assert len(C) == 3**t
    rows.append(
        {
            "length": n,
            "dimension": t,
            "catalog_types": len(codes),
            "generator_sha256": hashlib.sha256(
                codes.astype(np.uint8).tobytes()
            ).hexdigest(),
        }
    )
assert sum(row["catalog_types"] for row in rows) == 282
print(json.dumps(rows, indent=2))
print("PASS: all 282 catalog matrices have the stated rank and self-orthogonality.")
