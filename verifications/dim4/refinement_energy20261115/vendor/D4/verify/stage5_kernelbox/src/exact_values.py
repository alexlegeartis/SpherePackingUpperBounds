"""Lemma 1123: exact rationals of a_new, U_new, its strict ceiling, and the (26,5) next floor.
Same arithmetic as ledger.py with lam1 = 3/10^6, lam5 = 1/10^6. Run from any folder; inputs come from paths.py."""
from fractions import Fraction as Q
from decimal import Decimal, getcontext
from pathlib import Path
import json, sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

getcontext().prec = 40
acc = json.loads(paths.CENSUS.read_text())
a = Q(532, 1000); a2 = a * a
QN = lambda N: N * N * (Q(1, 8) - a2 / 4) - N * (1 - a2)
LAM = {1: Q(3, 10 ** 6), 2: Q(1, 10 ** 6), 3: Q(1, 10 ** 6), 4: Q(1, 10 ** 6), 5: Q(1, 10 ** 6)}
U_RAW = Q(acc["exact_raw_objective"]); GS = Q(acc["graph_saving"])
# The project shared budget file is not public. The same cross-check uses the stage-2 acceptance.
assert Q(json.loads(paths.ACCEPTANCE.read_text())["exact_source_objective"]) == U_RAW  # stage 2 and ledger use the same U_raw
floors = sorted((Q(g["new_floor"]) + LAM[g["first_radius"]] * max(QN(g["N"]), Q(0)), g["N"], g["first_radius"])
                for g in acc["profile"])
a_new, N0, j0 = floors[0]
nxt = floors[1]
assert (N0, j0) == (24, 1) and (nxt[1], nxt[2]) == (26, 5) and floors[0][0] < nxt[0]
assert a_new == Q(acc["new_beta"]) + Q(3, 10 ** 6) * QN(24) == Q(acc["new_beta"]) + Q(65799, 1562500000)
U_new = U_RAW - a_new / 8 - GS
ceil = Q(629521681653, 10 ** 12)
assert U_new <= ceil and ceil - U_new < Q(1, 10 ** 12)
dec = lambda x: str(Decimal(x.numerator) / Decimal(x.denominator))
out = dict(status="COMPUTED_EXACT",
           Q_24=str(QN(24)), Q_26=str(QN(26)),
           a_new=str(a_new), a_new_decimal=dec(a_new), a_new_rule="new_beta + (3/10^6) Q_24, group (24,1)",
           next_floor_26_5=str(nxt[0]), next_floor_26_5_decimal=dec(nxt[0]),
           next_floor_rule="floor(26,5) + (1/10^6) Q_26",
           gap_next_minus_a_new_decimal=dec(nxt[0] - a_new),
           U_new=str(U_new), U_new_decimal=dec(U_new),
           U_new_strict_ceiling="629521681653/1000000000000",
           old_strict_ceiling=acc["strict_upper_bound"],
           improvement_of_ceiling=str(Q(acc["strict_upper_bound"]) - ceil))
paths.OUT.mkdir(parents=True, exist_ok=True)
(paths.OUT / "exact_values.json").write_text(json.dumps(out, indent=1) + "\n")
print({k: (v if len(v) < 80 else v[:40] + "...") for k, v in out.items()})
