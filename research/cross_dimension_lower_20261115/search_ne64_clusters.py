"""Discovery: grow actual Ne64 dual normals and search multi-point clusters.

Finite pools, floating point enumeration and greedy cliques are discovery only.
Any saved improvement must pass the separate exact candidate auditor.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import argparse
import itertools
import json
import math
import time
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "research/novel_routes_20261113/dual64"


def points(H, bound, cap=500, node_cap=100000):
    R = np.linalg.cholesky(H).T
    z = np.zeros(len(H), dtype=np.int64)
    found = []
    nodes = 0
    capped = False
    def walk(i, left):
        nonlocal nodes, capped
        nodes += 1
        if nodes > node_cap or len(found) >= cap:
            capped = True
            return
        if i < 0:
            found.append(z.copy())
            return
        center = -R[i, i+1:]@z[i+1:]/R[i, i]
        radius = math.sqrt(max(0., left))/R[i, i]
        for x in range(math.ceil(center-radius-1e-8), math.floor(center+radius+1e-8)+1):
            z[i] = x
            walk(i-1, left-(R[i,i]*(x-center))**2)
            if capped:
                return
    walk(len(H)-1, bound+1e-7)
    return np.array(found, dtype=np.int64).reshape((-1, len(H))), capped


def cluster_score(Knum, rng, beta_cap=1.2):
    k = len(Knum)
    n = 64-k
    H = 3*np.linalg.inv(Knum)
    cov = 3**16*math.sqrt(np.linalg.det(Knum/3))
    best = ((3.)**(n/2)/cov, np.zeros((1,k),dtype=np.int64), 0.)
    W, capped = points(H, beta_cap)
    if capped:
        return best, True, len(W)
    q = np.einsum('vi,ij,vj->v', W, H, W)
    D = q[:,None]+q[None,:]-2*W@H@W.T
    zero = int(np.flatnonzero(np.all(W==0,axis=1))[0])
    for beta in np.unique(np.round(q,10)):
        if beta < 1e-8 or beta > beta_cap:
            continue
        ids = np.flatnonzero(q <= beta+1e-8)
        adj = []
        for i in ids:
            adj.append(sum(1<<j for j,t in enumerate(ids) if i != t and D[i,t] <= beta+1e-8))
        root = int(np.flatnonzero(ids==zero)[0])
        seeds = sorted(range(len(ids)), key=lambda i: -adj[i].bit_count())[:8]
        for seed in seeds:
            chosen = [root]
            opts = adj[root]
            if seed != root:
                chosen.append(seed)
                opts &= adj[seed]
            while opts:
                available = [j for j in range(len(ids)) if opts>>j&1]
                degrees = [(opts&adj[j]).bit_count() for j in available]
                maximum = max(degrees)
                tie = [j for j,d in zip(available,degrees) if d==maximum]
                j = int(rng.choice(tie))
                chosen.append(j)
                opts &= adj[j]
            labels = W[ids[chosen]]
            diameter = float(D[np.ix_(ids[chosen],ids[chosen])].max())
            score = len(labels)*((12-diameter)/4)**(n/2)/cov
            if score > best[0]:
                best = (score, labels, diameter)
    return best, False, len(W)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--limit", type=int, default=1000)
    ap.add_argument("--rank", type=int, default=6)
    ap.add_argument("--seed", type=int, default=581115)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    d = np.load(PARENT / "reduced_dual.npz")
    G,H,T,R = (d[key].astype(np.int64) for key in ("G","H","T","R"))
    V = np.load(PARENT / "short_pool.npz")["V"].astype(np.int64)
    V = np.vstack((V,-V))
    A5 = np.loadtxt(ROOT/"constructions/dim59/data/normal_coefficients.txt", dtype=np.int64)
    # Work directly in original dual coordinates, avoiding an isometry guess.
    A_pool = V@T
    seed_normals = [A5]
    best_pool = []
    report = {"rank":args.rank,"dimension":64-args.rank,"pool_signed":len(V),
              "seed":args.seed,"evaluated":0,"point_cap_skips":0,
              "improvements":[],"search_complete":False,
              "proof_scope":"Floating discovery only; finite pool and heuristic cliques."}
    rng = np.random.default_rng(args.seed)
    baseline = {58: 2**18/math.sqrt(3), 57:151706.083320,59:274522.5897411497}
    target = baseline.get(64-args.rank, float("inf"))
    seen = set()
    start = time.monotonic()
    for k in range(6,args.rank+1):
        target = baseline.get(64-k, float("inf"))
        candidates = []
        for A in seed_normals:
            K0 = A@H@A.T
            for v in A_pool:
                ip = A@H@v
                K = np.block([[K0,ip[:,None]],[ip[None,:],np.array([[12]])]])
                sig = K.tobytes()
                if sig in seen:
                    continue
                seen.add(sig)
                if np.linalg.eigvalsh(K.astype(float))[0] < 1e-7:
                    continue
                candidates.append((np.linalg.det(K),np.vstack((A,v)),K))
        candidates.sort(key=lambda item:item[0])
        candidates = candidates[:args.limit]
        best_pool = []
        for _,A,K in candidates:
            (score,S,beta),capped,count = cluster_score(K,rng)
            report["evaluated"] += 1
            report["point_cap_skips"] += int(capped)
            best_pool.append((score,A,K))
            if score > target*(1+1e-8):
                exactA = sp.Matrix(A.tolist())
                snf = smith_normal_form(exactA,domain=sp.ZZ)
                if any(abs(snf[i,i])!=1 for i in range(k)):
                    continue
                name=f"ne64_n{64-k}_s{len(S)}_{report['evaluated']}.npz"
                np.savez_compressed(args.output/name,G=G,H=H,A=A,Knum=K,S=S)
                row={"file":name,"dimension":64-k,"cosets":len(S),"density_float":score,
                     "beta_float":beta,"points_in_discovery_pool":count,"exact_audit_pending":True}
                report["improvements"].append(row)
                print("CANDIDATE",json.dumps(row),flush=True)
                target=score
            if report["evaluated"]%50==0:
                print("progress",k,report["evaluated"],"best",max(x[0] for x in best_pool),
                      "seconds",time.monotonic()-start,flush=True)
        best_pool.sort(key=lambda item:-item[0])
        seed_normals=[row[1] for row in best_pool[:6]]
        if best_pool:
            report[f"best_rank_{k}"]={"density_float":best_pool[0][0],"normal_gram_numerator":best_pool[0][2].tolist()}
    report["seconds"]=time.monotonic()-start
    (args.output/"search.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report),flush=True)


if __name__=="__main__":
    main()
