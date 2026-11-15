"""Exploratory normal-Gram optimization; no ambient embedding is presumed.

Optimize smaller subsets of the saved nine-point root clusters. Rational
minimum and embedding checks are still required for any resulting candidate.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import itertools
import json
import math
import argparse
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "novel_routes_20261113/higher_parent"


def short_vectors(K, bound):
    R = np.linalg.cholesky(K).T
    z = np.zeros(len(K), dtype=np.int64)
    out = []
    def walk(i, partial):
        if i < 0:
            if np.any(z): out.append(z.copy())
            return
        center = -R[i, i+1:] @ z[i+1:] / R[i, i]
        radius = math.sqrt(max(0., bound - partial)) / R[i, i]
        for x in range(math.ceil(center-radius-1e-8), math.floor(center+radius+1e-8)+1):
            z[i] = x
            term = (R[i, i] * (x-center))**2
            if partial+term <= bound+1e-7: walk(i-1, partial+term)
    walk(len(K)-1, 0.)
    return np.array(out, dtype=np.int64).reshape((-1, len(K)))


def optimize(k, outpath, circuits=False):
    if circuits:
        f = PARENT / f"gamma72_simplex_rank{k}.npz"
        K0 = np.load(f)["K"] if f.exists() else 8*np.eye(k,dtype=np.int64)
        S0 = np.vstack((np.zeros((1,k),dtype=np.int64),np.eye(k,dtype=np.int64)))
    else:
        a = np.load(PARENT / f"gamma72_root_cluster_rank{k}.npz")
        K0, S0 = a["K"], a["S"]
    pairs = list(itertools.combinations_with_replacement(range(k), 2))
    pcount = len(pairs)
    E = np.zeros((pcount, k, k))
    for t, (i, j) in enumerate(pairs): E[t, i, j] = E[t, j, i] = 1
    W0 = np.array(list(itertools.product(range(-1, 2), repeat=k)), dtype=np.int64)
    W0 = np.unique(np.vstack((W0[np.any(W0,axis=1)], short_vectors(K0, 8.01))), axis=0)
    n = 72-k
    rows = []
    if circuits:
        subsets = []
        for plus in range(1,k+1):
            for minus in range(k-plus+1):
                if plus+minus < 2: continue
                v = np.array([1]*plus+[-1]*minus+[0]*(k-plus-minus))
                subsets.append(([plus,minus],np.vstack((S0,v))))
    else:
        subsets = [(None, S0)] + [(i, np.delete(S0, i, axis=0)) for i in range(len(S0))]
    for removed, S in subsets:
        W = W0.copy()
        Z = np.array([x-y for i, x in enumerate(S) for y in S[:i]])
        def unpack(x): return np.einsum('t,tij->ij', x[:pcount], E)
        def inverse(x): return np.linalg.inv(unpack(x))
        def diameter(K):
            H = np.linalg.inv(K)
            return np.einsum('bi,ij,bj->b', Z, H, Z).max()
        def objective(x):
            sign, ld = np.linalg.slogdet(unpack(x))
            if sign <= 0 or x[-1] >= 8: return 1e8
            return .5*ld - n/2*math.log((8-x[-1])/4) - math.log(len(S))
        def jacobian(x):
            H = inverse(x)
            return np.r_[.5*np.einsum('ij,tij->t', H, E), n/2/(8-x[-1])]
        def cluster(x):
            H = inverse(x)
            return x[-1]-np.einsum('bi,ij,bj->b', Z, H, Z)
        def cluster_jac(x):
            A = Z @ inverse(x)
            return np.c_[np.einsum('bi,tij,bj->bt', A, E, A), np.ones(len(Z))]
        def eigen(x): return np.linalg.eigvalsh(unpack(x))[0]-.005
        def eigen_jac(x):
            v = np.linalg.eigh(unpack(x))[1][:,0]
            return np.r_[np.einsum('i,tij,j->t',v,E,v),0.]
        x = np.r_[[K0[i,j] for i,j in pairs], diameter(K0)]
        for cut in range(10):
            F = np.c_[np.einsum('bi,tij,bj->bt',W,E,W),np.zeros(len(W))]
            result = minimize(objective, x, jac=jacobian, method="SLSQP",
                bounds=[(8.,24.) if i==j else (-12.,12.) for i,j in pairs]+[(.01,2.)],
                constraints=[{"type":"ineq","fun":lambda x:F@x-8,"jac":lambda x:F},
                             {"type":"ineq","fun":cluster,"jac":cluster_jac},
                             {"type":"ineq","fun":eigen,"jac":eigen_jac}],
                options={"maxiter":600,"ftol":1e-10})
            x = result.x
            K = unpack(x)
            if np.linalg.eigvalsh(K)[0] <= 0: break
            new = short_vectors(K, 8-1e-5)
            if not len(new): break
            W = np.unique(np.vstack((W,new)),axis=0)
        row = {"rank":k,"removed_label":removed,"cosets":len(S),
               "optimizer_success":bool(result.success),"message":result.message,
               "density_relaxation":math.exp(-objective(x)),"beta":float(x[-1]),
               "normal_gram_float":K.tolist(),"labels":S.tolist(),
               "integral_embedded_construction":False,"cuts":cut,
               "minimum_vector_constraints":len(W)}
        rows.append(row)
        outpath.write_text(json.dumps(rows,indent=2)+"\n")
        print(k,removed,len(S),row["optimizer_success"],row["density_relaxation"],row["beta"],flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--rank",type=int,default=6)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--circuits",action="store_true")
    args = p.parse_args()
    optimize(args.rank,args.output,args.circuits)
