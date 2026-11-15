"""Lemma 1044 (kernelbox): fail-closed adaptive bisection proving
-K(y,z) >= lam * q(t), q = t^2 (t^2 - a^2), a = 532/1000, on
B_j = {t in [-1,-a], |y|,|z| in [1,R_j]}.

Deterministic. Root boxes: radius regions [R_k, R_{k+1}] (k < j) squared, times
t-chunks. Each root is bisected at midpoints until kcert.F_lower >= 0 on every
leaf; a box narrower than the floor widths fails the run. The split tree of
every root is saved so check_kernelbox.py can replay it without search.
Long-job control (project only): with --participant the parent polls
`board.py control` every 60 s and stops its pool on a nonzero exit. Run with nice -n 10.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, json, os, subprocess, sys, time
import multiprocessing as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kcert as kc
import kformula as kf

AA = Q(532, 1000)
MIN_W = (Q(1, 10 ** 7), Q(1, 10 ** 7), Q(1, 10 ** 7))


def roots(j, nt):
    regs = [(kf.BREAKS[k], kf.BREAKS[k + 1], k > 0) for k in range(j)]
    ts = [Q(-1) + (Q(-1) * -1 - AA) * Q(i, nt) for i in range(nt + 1)]
    out = []
    for ra in regs:
        for rb in regs:
            for i in range(nt):
                out.append(((ra[0], ra[1]), (rb[0], rb[1]), (ts[i], ts[i + 1]), ra[2], rb[2]))
    return out


def solve(args):
    root, lam = args
    br, bs, bt, open_r, open_s = root
    stats = dict(leaves=0, evals=0, min_lower=None, fail=None)

    def rec(box, depth):
        f = kc.F_box(box, lam, AA, open_r, open_s)
        stats["evals"] += 1
        lo = f.lower()
        if lo >= 0:
            stats["leaves"] += 1
            v = Q(*lo.mid().man_exp()) if False else float(lo.mid())
            stats["min_lower"] = v if stats["min_lower"] is None else min(stats["min_lower"], v)
            return "L"
        sc = kc.split_scores(f)
        order = sorted(range(3), key=lambda i: -sc[i])
        for d in order:
            a, b = box[d]
            if b - a > MIN_W[d] * 2:
                break
        else:
            stats["fail"] = [[str(x) for x in iv] for iv in box]
            return "F"
        m = (a + b) / 2
        left = list(box); right = list(box)
        left[d] = (a, m); right[d] = (m, b)
        L = rec(tuple(left), depth + 1)
        if stats["fail"]:
            return [d, L, "?"]
        R = rec(tuple(right), depth + 1)
        return [d, L, R]

    tree = rec((br, bs, bt), 0)
    enc = dict(r=[str(x) for x in br], s=[str(x) for x in bs], t=[str(x) for x in bt],
               open_r=open_r, open_s=open_s, tree=tree)
    return enc, stats


def control_ok(pid):
    if not pid:  # public package: no project board, no remote stop control
        return True
    res = subprocess.run([sys.executable, "board/board.py", "control", "--participant", pid],
                         capture_output=True, text=True)
    return res.returncode == 0


def atomic(path, obj):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj))
    os.replace(tmp, path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--j", type=int, required=True)
    ap.add_argument("--lam", required=True)
    ap.add_argument("--nt", type=int, default=94)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", required=True)
    ap.add_argument("--participant", default="", help="project board stop control; empty disables it")
    ap.add_argument("--only", type=int, nargs="*", help="debug: indices of root boxes")
    args = ap.parse_args()
    lam = Q(args.lam)
    rs = roots(args.j, args.nt)
    if args.only:
        rs = [rs[i] for i in args.only]
    out = Path(args.out); status = out.with_suffix(".status.json")
    atomic(status, dict(state="running", pid=os.getpid(), roots=len(rs), started=time.time()))
    t0 = time.time(); results = [None] * len(rs); stopped = False
    with mp.Pool(args.workers) as pool:
        handles = [pool.apply_async(solve, ((r, lam),)) for r in rs]
        last = time.time()
        while not all(h.ready() for h in handles):
            time.sleep(2)
            if time.time() - last > 60:
                last = time.time()
                done = sum(h.ready() for h in handles)
                atomic(status, dict(state="running", pid=os.getpid(), roots=len(rs), done=done,
                                    elapsed=time.time() - t0))
                if not control_ok(args.participant):
                    stopped = True; pool.terminate(); break
        if not stopped:
            for i, h in enumerate(handles):
                results[i] = h.get()
    if stopped:
        atomic(status, dict(state="stopped_by_control", elapsed=time.time() - t0))
        sys.exit(3)
    fails = [s["fail"] for _, s in results if s["fail"]]
    summary = dict(j=args.j, lam=str(lam), a=str(AA), nt=args.nt, roots=len(rs),
                   leaves=sum(s["leaves"] for _, s in results), evals=sum(s["evals"] for _, s in results),
                   min_leaf_lower=min(s["min_lower"] for _, s in results if s["min_lower"] is not None),
                   failures=fails, elapsed=time.time() - t0,
                   status="ALL_LEAVES_NONNEGATIVE" if not fails else "FAILED")
    out.write_text(json.dumps(dict(summary=summary, roots=[e for e, _ in results])))
    atomic(status, dict(state="done", **{k: v for k, v in summary.items() if k != "failures"}, nfail=len(fails)))
    print(json.dumps({k: v for k, v in summary.items() if k != "failures"}), "nfail", len(fails))
    sys.exit(0 if not fails else 1)
