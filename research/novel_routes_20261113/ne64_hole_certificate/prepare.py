"""Prepare exact outward enclosures for a proposed deep Ne64 dual coset.

This establishes the exact input identities. It does not by itself certify
the hole: every shard of interval_enum must complete without a witness.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp

P = Path(__file__).resolve().parent
SOURCE = P.parent / "higher_parent"
CLASS_ID = 4
SCALE = 2**40


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    data = np.load(SOURCE / "ne64_primal_best.npz")
    G, H, T, Ti, Q = (sp.Matrix(data[k].tolist())
                       for k in ("G", "H", "T", "Ti", "R"))
    base = np.load(P.parent / "dual64" / "source.npz")
    assert G == sp.Matrix(base["G"].tolist())
    assert H == sp.Matrix(base["H"].tolist())
    assert G * H == 3 * sp.eye(64)
    assert T * Ti == sp.eye(64)
    assert Q == T * G * T.T
    assert abs(T.det()) == 1
    assert G.det() == 3**32
    assert G == G.T and all(G[i, i] % 2 == 0 for i in range(64))

    candidate = next(x for x in json.loads((SOURCE / "coset_classes.json").read_text())
                     if x["id"] == CLASS_ID)
    a = sp.Matrix([candidate["original_residue"]])
    c = sp.Matrix([candidate["reduced_residue"]])
    assert all(int(x) in (0, 1, 2) for x in a)
    assert all(int(x) in (0, 1, 2) for x in c)
    assert all(int(x) % 3 == 0 for x in c * T - a)
    assert all(int(x) % 3 == 0 for x in a * G)
    assert any(x for x in a)
    assert int((a * G * a.T)[0]) % 18 == 6
    # Therefore h=a/3 lies in L*, has order exactly 3 modulo L, and
    # every representative of h+L has norm congruent to 2/3 modulo 2.
    # Absence of integer z == c (mod 3) with z Q z^t <= 60 proves
    # min(h+L) >= 78/9 = 26/3.

    L, D = Q.LDLdecomposition(hermitian=False)
    assert L * D * L.T == Q
    assert all(D[i, i] > 0 for i in range(64))
    with (P / "ldl_bounds.txt").open("w") as f:
        for i in range(64):
            for x in [D[i, i]] + [L[j, i] for j in range(64)]:
                lo, hi = int(sp.floor(x * SCALE)), int(sp.ceiling(x * SCALE))
                assert abs(lo) < 2**53 and abs(hi) < 2**53
                assert sp.Rational(lo, SCALE) <= x <= sp.Rational(hi, SCALE)
                f.write(f"{lo} {hi}\n")
    np.savetxt(P / "gram.txt", np.array(Q.tolist(), dtype=np.int64), fmt="%d")
    (P / "coset.txt").write_text(" ".join(map(str, candidate["reduced_residue"])) + "\n")
    setup = {
        "dimension": 64,
        "class_id": CLASS_ID,
        "scale": SCALE,
        "inclusive_integer_norm_bound": 60,
        "congruence_modulus": 3,
        "residue_norm_mod_18": 6,
        "desired_coset_minimum": "26/3",
        "exact_basis_change_verified": True,
        "exact_ldl_identity_verified": True,
        "exact_dual_and_coset_type_verified": True,
        "parent_minimum_12": "inherited from published Ne64 theorem; not enumerated here",
        "original_residue": candidate["original_residue"],
        "reduced_residue": candidate["reduced_residue"],
        "sha256": {name: digest(P / name) for name in ("ldl_bounds.txt", "gram.txt", "coset.txt")},
        "source_sha256": digest(SOURCE / "ne64_primal_best.npz"),
        "hole_certified_by_setup_alone": False,
    }
    (P / "setup.json").write_text(json.dumps(setup, indent=2) + "\n")
    print(json.dumps({k: v for k, v in setup.items()
                      if k not in ("original_residue", "reduced_residue", "sha256")}, indent=2))


if __name__ == "__main__":
    main()
