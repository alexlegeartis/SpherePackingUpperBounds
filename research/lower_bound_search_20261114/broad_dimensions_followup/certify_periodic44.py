"""Verified positive-stress existence bound for the sampled 44D contacts.

Binary64 dot-product error bounds are widened with epsilon=2^-52 (twice the
unit roundoff). The proof uses exact integer contact equations, not rounded
orthonormal coordinates. See README for the Neumann-series argument.
"""

from pathlib import Path
import json, time
import numpy as np
from flint import fmpz_mat

P = Path(__file__).resolve().parent
started = time.time()
z = np.load(P / "antipode44_0.npz")
d = np.load(P / "periodic44_stress_inputs.npz")
W = d["W"]
labels = d["labels"]
weights = d["weights"].copy()
Q = z["Q"]
ii, jj = np.triu_indices(44)
p = len(ii)
C = np.zeros((p + 176, len(W)), dtype=np.int64)
C[:p] = (W[:, ii] * W[:, jj]).T
# Reconstruct the primitive section and every contact without trusting the
# random root sampler or the LP. The ambient minimum six is a literature input.
G = np.load(P / "ambient48.npz")["G"]
F = fmpz_mat(G.tolist())
assert F.det() == 1
Gi = np.array(F.inv().tolist(), dtype=np.int64)
assert np.array_equal(G @ Gi, np.eye(48, dtype=np.int64))
assert abs(fmpz_mat(z["U"].tolist()).det()) == 1
assert np.array_equal(Q, z["U"] @ Gi[np.ix_(z["J"], z["J"])] @ z["U"].T)
assert fmpz_mat(Q.tolist()).det() == 605
assert all(fmpz_mat(Q[:k, :k].tolist()).det() > 0 for k in range(1, 45))
K = z["signs"][:, None] * G[np.ix_(z["indices"], z["indices"])] * z["signs"][None, :]
assert np.array_equal(K, z["K"])
assert np.array_equal(K @ z["M55"], 55 * np.eye(4, dtype=np.int64))
assert np.array_equal(z["U"] @ z["Uinv"], np.eye(44, dtype=np.int64))
Csection = z["signs"][:, None] * G[np.ix_(z["indices"], z["J"])]
assert np.array_equal(Csection, z["C"])
assert np.array_equal(
    z["shifts55"],
    np.r_[np.zeros((1, 44), dtype=np.int64), -z["M55"] @ Csection] @ z["Uinv"],
)
assert np.all((W - z["shifts55"][labels[:, 1]] + z["shifts55"][labels[:, 0]]) % 55 == 0)
assert 44**2 * int(abs(Q).max()) * int(abs(W).max()) ** 2 < 2**63
assert np.all(np.einsum("vi,ij,vj->v", W, Q, W) == 17270)
assert set(map(tuple, labels)) <= set(map(tuple, z["active_pairs"]))
for a in range(1, 5):
    C[p + (a - 1) * 44 : p + a * 44] = (
        W * ((labels[:, 1] == a).astype(int) - (labels[:, 0] == a).astype(int))[:, None]
    ).T
inv = fmpz_mat(Q.tolist()).inv()
targetq = [[x * 3025 for x in row] for row in inv.tolist()]
assert all(x.denominator == 1 for row in targetq for x in row)
target = np.r_[np.array(targetq, dtype=np.int64)[ii, jj], np.zeros(176, dtype=np.int64)]
# Every integer dot product has sum of absolute values < 2^53, by Cauchy-
# Schwarz and these exact row square sums. Every multiply and add in an IEEE
# binary64 dot product is therefore exact, in any parenthesization.
max_entry = int(np.max(abs(C)))
assert len(W) * max_entry * max_entry < 2**63
diagonal = np.sum(C * C, axis=1, dtype=np.int64)
assert int(diagonal.max()) < 2**53
Cd = C.astype(float)
G = Cd @ Cd.T
assert np.array_equal(G, G.T) and np.array_equal(np.diag(G), diagonal.astype(float))
assert np.array_equal(G, np.rint(G))
powers = np.ceil(np.log2(np.maximum(np.max(abs(C), axis=1), 1))).astype(int)
scale = np.ldexp(np.ones(len(C)), -powers)
A = Cd * scale[:, None]
t = target * scale
Gs = G * scale[:, None] * scale[None, :]
H = np.linalg.inv(Gs)
print(
    "formed exact integer Gram",
    len(C),
    len(W),
    "max square sum",
    int(diagonal.max()),
    flush=True,
)
for iteration in range(3):
    weights = weights + A.T @ (H @ (t - A @ weights))
eps = 2.0**-52


def gamma(k):
    return (k * eps) / (1 - k * eps)


def inflate(x):
    return np.nextafter(x * (1 + 1e-11) + 1e-280, np.inf)


gd = gamma(len(C))
gn = gamma(len(W))
I = np.eye(len(C))
B = Gs @ H
absproduct = abs(Gs) @ abs(H)
Ebound = abs(I - B) + gd / (1 - gd) * absproduct + eps * (I + abs(B))
eta = float(inflate(np.max(np.sum(Ebound, axis=1)) / (1 - gd)))
R = A.T @ H
Rabs = abs(A.T) @ abs(H)
right_inverse_norm = float(
    inflate(np.max(np.sum(abs(R) + gd / (1 - gd) * Rabs, axis=1)) / (1 - gd))
)
value = A @ weights
absolute = abs(A) @ abs(weights)
residual = float(
    inflate(
        np.max(abs(t - value) + gn / (1 - gn) * absolute + eps * (abs(t) + abs(value)))
    )
)
correction = (
    float(inflate(right_inverse_norm * residual / (1 - eta)))
    if eta < 1
    else float("inf")
)
margin = float(weights.min()) - correction
result = {
    "contacts": len(W),
    "equations": len(C),
    "max_exact_Gram_diagonal": str(int(diagonal.max())),
    "neumann_residual_bound": eta,
    "right_inverse_norm_bound": right_inverse_norm,
    "stress_residual_bound": residual,
    "weight_min": float(weights.min()),
    "weight_max": float(weights.max()),
    "correction_bound": correction,
    "positive_margin": margin,
    "verified": bool(eta < 1 and margin > 0),
    "seconds": time.time() - started,
}
print(json.dumps(result), flush=True)
assert result["verified"]
