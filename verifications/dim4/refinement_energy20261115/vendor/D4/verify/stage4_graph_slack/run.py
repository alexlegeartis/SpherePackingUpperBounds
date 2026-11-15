"""Stage 4: broad antipodal graph capacity 1208 and the density saving 3/8000000. Fail closed.

Runs, each as a child process (nonzero exit, timeout or missing output is FAIL):
  1 src/verify_annular_count_151.py  exact degree-20 polynomial and Bernstein sign
    proof: at most 151 neighbours with radii in [33/28, 7/4].
  2 src/verify_projected_eight.py  exact Legendre polynomial and Bernstein sign
    proof: at most 8 projected neighbours on the full interval [-1, 223/735].
  3 src/gap_check.py  Q(7/8,7/8) = 223/735, 7/8 > 41/49, 151*8 = 1208,
    (1358-1208)/(8*50000000) = 3/8000000.
Then this script rechecks every exact step itself, including the corner
cosines and the place of the graph charge in U_raw, from data/v7/origin.json,
data/v7/graph_energies.json, the census floors and the stage-2 acceptance.

Usage: python3 run.py --out DIR [--census DIR] [--stage2 DIR]
--census / --stage2 take fresh upstream output folders; their RESULT.json must say PASS.
Writes DIR/RESULT.json; exit 0 iff status PASS. Standard library only.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time

sys.set_int_max_str_digits(0)
STAGE = Path(__file__).resolve().parent
SRC = STAGE / "src"
PKG = STAGE.parents[1]
DATA = PKG / "data"
TIMEOUT = 600


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def child(args, log, name):
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable] + args, capture_output=True, text=True, timeout=TIMEOUT, cwd=STAGE,
                           env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"{name}: timeout")
    log[name] = dict(exit=r.returncode, seconds=round(time.time() - t0, 2), stdout_tail=r.stdout[-800:], stderr_tail=r.stderr[-400:])
    if r.returncode != 0:
        raise RuntimeError(f"{name}: exit {r.returncode}: {(r.stdout + r.stderr)[-800:]}")
    return r.stdout


def upstream(census, stage2):
    """Resolve the census and stage-2 files. A folder given on the command line is a fresh
    upstream output: its RESULT.json must say PASS for the right stage, or this stage fails."""
    def passed(d, stage):
        r = d / "RESULT.json"
        if not r.exists():
            raise RuntimeError(f"missing {r}")
        res = json.loads(r.read_text())
        if res.get("stage") != stage or res.get("status") != "PASS":
            raise RuntimeError(f"{r}: stage {res.get('stage')!r} status {res.get('status')!r}, need {stage!r} PASS")
    if census:
        d = Path(census).resolve(); passed(d, "stage3_census")
        names = [n for n in ("census_floors.json", "FINITE_COUNT_ACCEPTANCE.json") if (d / n).exists()]
        if not names:
            raise RuntimeError(f"no census_floors.json or FINITE_COUNT_ACCEPTANCE.json in {d}")
        census_file = d / names[0]
    else:
        census_file = DATA / "census_floors.json"
    if stage2:
        d = Path(stage2).resolve(); passed(d, "stage2_origin")
    else:
        d = PKG / "verify/stage2_origin/reference"
    return census_file, d / "ACCEPTANCE.json"


def h(r, s):
    """Largest cosine between unit-separated points at radii r, s: (r^2+s^2-1)/(2rs)."""
    return (r * r + s * s - 1) / (2 * r * s)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--census", help="fresh stage 3 output folder: RESULT.json PASS plus census_floors.json or FINITE_COUNT_ACCEPTANCE.json (default: data/census_floors.json)")
    ap.add_argument("--stage2", help="fresh stage 2 output folder: RESULT.json PASS plus ACCEPTANCE.json (default: verify/stage2_origin/reference)")
    args = ap.parse_args()
    out = Path(args.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    res = dict(stage="stage4_graph_slack", status="FAIL"); log = {}; t0 = time.time()
    try:
        census_file, acc_file = upstream(args.census, args.stage2)
        res["census_file"] = str(census_file); res["stage2_file"] = str(acc_file)
        cert = SRC / "annular_count_151_certificate.json"
        projected = DATA / "projected_eight.json"
        inputs = [cert, projected, DATA / "v7/origin.json", DATA / "v7/graph_energies.json", census_file, acc_file]
        res["inputs"] = {str(p.relative_to(PKG)) if p.is_relative_to(PKG) else str(p): sha(p) for p in inputs}

        # 1 the 151-neighbour count, recomputed from the polynomial certificate
        a151 = json.loads(child([str(SRC / "verify_annular_count_151.py"), str(cert)], log, "annular_count_151"))
        assert a151["status"] == "EXACT_ANNULAR_COUNT_151_PASS" and a151["count_bound"] == 151
        assert Q(a151["exact_count_upper"]) < 152 and a151["spherical_cutoff"] == "41/49"
        sp = a151["whole_interval_sign_proof"]
        (out / "annular_count_151.json").write_text(json.dumps(a151, indent=1) + "\n")

        # 2 the exact projected S^2 code certificate
        projected_result = json.loads(child(
            [str(SRC / "verify_projected_eight.py"), str(projected)],
            log, "projected_eight"))
        assert projected_result["status"] == "EXACT_PROJECTED_EIGHT_PASS"
        assert projected_result["count_bound"] == 8
        assert Q(projected_result["cutoff"]) == Q(223, 735)
        assert Q(projected_result["objective"]) < 9
        (out / "projected_eight.json").write_text(json.dumps(projected_result, indent=1) + "\n")

        # 3 the exact geometric gap script
        gap = child([str(SRC / "gap_check.py")], log, "gap_check").split()
        assert gap == ["Q(7/8,7/8)=223/735", "gap=22/735", "7/8-41/49=15/392", "ordered=1208", "conditional_saving=3/8000000"], gap

        # 4 every exact geometric step again, here
        L, R = Q(33, 28), Q(7, 4)
        corners = dict(LL=h(L, L), LR=h(L, R), RR=h(R, R))
        assert corners == dict(LL=Q(697, 1089), LR=Q(41, 49), RR=Q(41, 49))
        rho, beta = max(corners.values()), Q(7, 8)
        assert R - 1 / R == L                                   # 33/28 = R - 1/R keeps the mixed corner at 41/49
        Qbb = (rho - beta * beta) / (1 - beta * beta)           # largest projected cosine
        assert Qbb == Q(223, 735) and Q(1, 3) - Qbb == Q(22, 735) > 0
        assert rho - beta < 0                                   # dQ/dh_i has numerator rho h_i - h_j <= rho - beta < 0
        assert beta > rho                                       # zero-projection vertex has degree <= 1
        degree_cap = projected_result["count_bound"]
        edges = 151 * degree_cap
        assert edges == 1208
        # selected edges are exterior: |y-z|^2 >= 2 L^2 (1 + 7/8) > (7/4)^2
        assert 2 * L * L * (1 + beta) == Q(16335, 3136) > Q(49, 16)

        # 5 where the graph charge sits in the objective
        origin = json.loads((DATA / "v7/origin.json").read_text())
        graphs = json.loads((DATA / "v7/graph_energies.json").read_text())
        w = list(map(Q, graphs["weights"]))
        assert w[:2] == [0, 0] and w[2] == Q(1, 50000000)
        charges = [(w[2], 1358), (Q(graphs["outer_graph_weight"]), 120), (Q(graphs["narrow_graph_weight"]), 604),
                   (Q(graphs["moderate_graph_weight"]), 1880)]
        charge = sum(x * c for x, c in charges)
        assert charge == Q(4769, 25000000)
        U_raw = Q(origin["objective"]) + charge / 8
        census = json.loads(census_file.read_text()); acc = json.loads(acc_file.read_text())
        assert U_raw == Q(census["exact_raw_objective"]) == Q(acc["exact_source_objective"])
        saving = w[2] * (1358 - edges) / 8
        assert saving == Q(3, 8000000) == Q(census["graph_saving"])
        assert Q(census["exact_objective"]) == U_raw - Q(census["new_beta"]) / 8 - saving
        res.update(status="PASS", annular_count=151, exact_count_upper=a151["exact_count_upper"],
                   projected_eight=projected_result,
                   annular_polynomial_degree=20, bernstein_nodes=sp["nodes"], bernstein_leaves=sp["leaves"],
                   bernstein_depth=sp["depth"], corner_cosines={k: str(v) for k, v in corners.items()},
                   rho=str(rho), beta=str(beta), projected_cosine_max=str(Qbb), gap_to_one_third=str(Q(1, 3) - Qbb),
                   degree_cap=degree_cap, ordered_edge_capacity=edges, old_capacity=1358,
                   broad_graph_weight=str(w[2]), graph_origin_charge=str(charge), U_raw=str(U_raw),
                   slack_in_energy=str(w[2] * (1358 - edges)), density_saving=str(saving),
                   classical_nine_point_theorem="historical alternative; not used by this replay")
    except Exception as e:  # fail closed
        res["status"] = "FAIL"; res["error"] = repr(e)[:3000]
    res["children"] = log
    res["seconds"] = round(time.time() - t0, 2)
    (out / "RESULT.json").write_text(json.dumps(res, indent=1) + "\n")
    print(res["status"], res.get("error", ""))
    sys.exit(0 if res["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
