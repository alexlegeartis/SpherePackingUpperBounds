"""Lemma 1099 (hrlemma, D4-12 ledger route): exact checks for the count-only far-pair slack lemma.
Run from any folder:
    python3 src/hrlemma_check.py
Exit 0 iff every check passes. Writes paths.OUT/hrlemma_check.json. Exact rationals only.
Checks: (1) Gegenbauer (parameter 1, U_k) expansion of q and the signs of its coefficients;
(2) cap s_j = 1 - 1/(2 R_j^2) <= a and the corner inequality behind the cap bound, j = 1..5;
(3) far pairs: 2 + 2a > 49/16; (4) Q_N table and positivity threshold;
(5) the cited kernel certificates are present with status ALL_LEAVES_NONNEGATIVE (read only,
not re-verified here; check_kernelbox.py gate 4 replays every leaf).
"""
from fractions import Fraction as Q
from pathlib import Path
import json, sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths
ROOT = paths.PKG
HERE = paths.OUT
KB = paths.KBOX
kernel = json.loads((paths.V7 / "kernel.json").read_text())
R = [Q(x) for x in kernel["breaks"][1:6]]                      # R_1..R_5
assert R == [Q(501473, 500000), Q(1005643, 1000000), Q(50379, 50000), Q(50767, 50000), Q(102597, 100000)]
a = Q(133, 250); a2 = a * a
out = {}; ok = True
def check(name, cond, value=None):
    global ok
    ok &= bool(cond); out[name] = dict(ok=bool(cond), value=value)
    print(("PASS " if cond else "FAIL ") + name, value if value is not None else "")

def U(k):
    p, c = [Q(1)], [Q(0), Q(2)]
    if k == 0: return p
    for _ in range(k - 1):
        n = [Q(0)] + [2 * x for x in c]
        for i, x in enumerate(p): n[i] -= x
        p, c = c, n
    return c
def padd(*ps):
    n = max(len(p) for p in ps); return [sum(p[i] if i < len(p) else 0 for p in ps) for i in range(n)]
def sc(p, x): return [x * y for y in p]
def ev(p, t): return sum(c * t ** i for i, c in enumerate(p))

# (1) q = t^4 - a^2 t^2 = c4 U4 + c2 U2 + c0 U0
q = [Q(0), Q(0), -a2, Q(0), Q(1)]
c4, c2, c0 = Q(1, 16), Q(3, 16) - a2 / 4, Q(1, 8) - a2 / 4
check("gegenbauer_identity", padd(sc(U(4), c4), sc(U(2), c2), sc(U(0), c0)) == q, [str(c4), str(c2), str(c0)])
check("coefficients_nonnegative", c4 >= 0 and c2 >= 0 and c0 > 0, [float(c4), float(c2), float(c0)])
check("q_at_1", ev(q, Q(1)) == 1 - a2, str(1 - a2))
check("c0_exact", c0 == Q(13561, 250000), str(c0))
# q <= 0 on [-a, a]: q = t^2 (t^2 - a^2)
check("q_nonpositive_on_[-a,a]", all(ev(q, t) <= 0 for t in (-a, Q(-1, 2), Q(0), Q(1, 2), a)) and True, "t^2 (t^2-a^2), factor sign")

# (2) caps
caps = []
for j, Rj in enumerate(R, 1):
    s = 1 - 1 / (2 * Rj * Rj)
    corner = (Rj - 1) * (Rj * Rj - Rj - 1) <= 0           # R/2 <= 1 - 1/(2R^2)
    caps.append(dict(j=j, R=str(Rj), s=str(s), s_decimal=float(s), s_le_a=s <= a, corner_ok=corner))
    check(f"cap_j{j}", s <= a and corner and Q(1, 2) <= s and 1 <= Rj, float(s))
out["caps"] = caps

# (3) far pairs
check("far_pairs", 2 + 2 * a > Q(49, 16), [str(2 + 2 * a), str(Q(49, 16))])

# (4) Q_N
def QN(N): return N * N * c0 - N * (1 - a2)
check("QN_formula", all(QN(N) == Q(N * (13561 * N - 179244), 250000) for N in range(0, 40)), "N(13561N-179244)/250000")
thr = min(N for N in range(1, 40) if QN(N) > 0)
check("QN_positive_from", thr == 14 and QN(13) < 0, thr)
table = {N: dict(exact=str(QN(N)), decimal=float(QN(N))) for N in range(14, 27)}
out["QN"] = table
for N in (24, 25, 26):
    print("   Q_%d = %s = %.6f" % (N, QN(N), float(QN(N))))
check("Q24_Q25_Q26", (QN(24), QN(25), QN(26)) == (Q(87732, 6250), Q(159781, 10000), Q(2253446, 125000)),
      [str(QN(24)), str(QN(25)), str(QN(26))])

# (5) cited kernel certificates (read only)
cert = {}
for j, lam in ((1, "3/1000000"), (5, "1/1000000")):
    # The public package has no run receipts (kbox_j*.status.json); the same fields are read
    # from the summary block of the published split tree.
    p = KB / f"kbox_j{j}.json"
    st = json.loads(p.read_text())["summary"] if p.exists() else {}
    good = st.get("status") == "ALL_LEAVES_NONNEGATIVE" and st.get("lam") == lam \
        and st.get("a") == "133/250" and st.get("j") == j and st.get("failures") == []
    cert[j] = dict(path=str(p.relative_to(ROOT)), status=st.get("status"), lam=st.get("lam"), leaves=st.get("leaves"))
    check(f"cited_kbox_j{j}", good, cert[j])
out["cited_certificates"] = cert
out["lambda"] = {"j=1": "3/1000000 (B_1)", "j=2..5": "1/1000000 (B_5 contains B_j)"}
HERE.mkdir(parents=True, exist_ok=True)
(HERE / "hrlemma_check.json").write_text(json.dumps(dict(status="PASS" if ok else "FAIL", checks=out), indent=1, default=str) + "\n")
print("ALL PASS" if ok else "FAILED")
sys.exit(0 if ok else 1)
