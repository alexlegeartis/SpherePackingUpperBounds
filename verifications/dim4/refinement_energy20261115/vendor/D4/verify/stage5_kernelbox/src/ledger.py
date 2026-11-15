"""Lemma 1044 (kernelbox) part (2): new universal floor and strict ceiling.

a_new = min over the 124 accepted groups g = (N, j) of floor_g + lam(j) * max(Q_N, 0),
lam(1) = LAM1 (certificate on B_1), lam(j) = LAM5 for j = 2..5 (B_j is inside B_5).
Q_N = N^2 (1/8 - a^2/4) - N (1 - a^2), a = 532/1000. H_R >= 0 always, hence max(., 0).
U_new = U_raw - a_new/8 - 3/8000000, rounded up to 12 decimals (strict ceiling).
Exact rationals throughout. Run from any folder; the census file comes from paths.py.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, json, sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

ap = argparse.ArgumentParser()
ap.add_argument("--lam1", required=True)
ap.add_argument("--lam5", default="0")
ap.add_argument("--out")
args = ap.parse_args()
LAM1, LAM5 = Q(args.lam1), Q(args.lam5)
acc = json.loads(paths.CENSUS.read_text())
a = Q(532, 1000); a2 = a * a
QN = lambda N: N * N * (Q(1, 8) - a2 / 4) - N * (1 - a2)
U_RAW = Q(acc["exact_raw_objective"]); GS = Q(acc["graph_saving"]); A_PLUS = Q(acc["new_beta"])
assert Q(acc["exact_objective"]) == U_RAW - A_PLUS / 8 - GS
assert GS == Q(3, 8000000) and len(acc["profile"]) == 124
assert min(Q(g["new_floor"]) for g in acc["profile"]) == A_PLUS
rows = []
for g in acc["profile"]:
    j = g["first_radius"]; assert j in (1, 2, 3, 4, 5)
    lam = LAM1 if j == 1 else LAM5
    inc = lam * max(QN(g["N"]), Q(0))
    rows.append((Q(g["new_floor"]) + inc, g["N"], j, inc))
rows.sort()
a_new = rows[0][0]
U_new = U_RAW - a_new / 8 - GS
D = 10 ** 12
ceil = Q(-((-U_new.numerator * D) // U_new.denominator), D)
assert ceil >= U_new and ceil - U_new < Q(1, D)
old = Q(acc["strict_upper_bound"])
res = dict(status="COMPUTED_EXACT_LEDGER", lam1=str(LAM1), lam5=str(LAM5), a=str(a),
           a_plus=float(A_PLUS), a_new=str(a_new), a_new_decimal=float(a_new), a_gain=float(a_new - A_PLUS),
           binding=[dict(N=N, j=j, new_floor=float(v), increment=float(i)) for v, N, j, i in rows[:6]],
           U_raw_decimal=float(U_RAW), U_new_exact_decimal=f"{float(U_new):.15f}",
           U_new_strict_ceiling=str(ceil), U_new_ceiling_decimal=f"{float(ceil):.12f}",
           old_strict_ceiling=str(old), improvement=str(old - ceil), improvement_decimal=float(old - ceil))
print(json.dumps(res, indent=1))
if args.out:
    Path(args.out).write_text(json.dumps(res, indent=1) + "\n")
