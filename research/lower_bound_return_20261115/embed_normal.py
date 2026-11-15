"""Discover an explicit Gamma72 embedding; verify every saved identity.

The finite pool and backtracking cap restrict discovery only. A found embedding
is checked by exact integer products and Smith invariant factors.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import argparse
import json
import time
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "novel_routes_20261113/higher_parent"


def pool(path, size, seed):
    G = np.load(PARENT / "G72.npy").astype(np.int64)
    if path.exists():
        a = np.load(path)
        assert np.array_equal(a["G"],G)
        return G,a["V"]
    E = np.eye(72,dtype=np.int16)
    V = [E]
    for f in sorted(PARENT.glob("gamma72_*rank*.npz")):
        a = np.load(f)
        if "B" in a.files: V.append(a["B"])
    for i in range(72):
        for j in range(i):
            if abs(G[i,j]) == 4: V.append((E[i]-int(np.sign(G[i,j]))*E[j])[None,:])
    V = np.vstack(V).astype(np.int16)
    norms = np.einsum('vi,ij,vj->v',V.astype(np.int64),G,V.astype(np.int64))
    V = V[norms==8]
    V *= np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)])[:,None]
    V = np.unique(V,axis=0)
    rng = np.random.default_rng(seed)
    for step in range(30):
        selected = rng.choice(len(V),min(350,len(V)),replace=False)
        C = V[selected].astype(np.int64) @ G @ V.T
        i,j = np.where(np.abs(C)==4)
        W = V[selected[i]]-np.sign(C[i,j])[:,None]*V[j]
        W = W.astype(np.int16)
        W *= np.sign(W[np.arange(len(W)),np.argmax(W!=0,axis=1)])[:,None]
        V = np.unique(np.vstack((V,W)),axis=0)
        if len(V)>size: V=V[rng.choice(len(V),size,replace=False)]
        print("pool",step,len(V),flush=True)
        if len(V)>=size: break
    assert np.all(np.einsum('vi,ij,vj->v',V.astype(np.int64),G,V.astype(np.int64))==8)
    np.savez_compressed(path,G=G,V=V)
    return G,V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--pool",type=Path,required=True)
    ap.add_argument("--size",type=int,default=100000)
    ap.add_argument("--seed",type=int,default=661115)
    ap.add_argument("--cap",type=int,default=30000)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    G,V=pool(args.pool,args.size,args.seed)
    V=np.vstack((V,-V)).astype(np.int64)
    VG=V@G
    cache={}
    def ips(i):
        if i not in cache:
            cache[i]=VG@V[i]
            if len(cache)>256:cache.pop(next(iter(cache)))
        return cache[i]
    candidates=json.loads(args.candidates.read_text())
    for ci,row in enumerate(candidates):
        K=np.array(row["K"],dtype=np.int64)
        # Strongest first constraints reduce the branching factor.
        order=[int(np.argmax(np.abs(K-np.diag(K.diagonal())).sum(axis=1)))]
        while len(order)<len(K):
            left=[i for i in range(len(K))if i not in order]
            order.append(max(left,key=lambda i:sum(abs(int(K[i,j])) for j in order)))
        Kp=K[np.ix_(order,order)]
        calls=0;start=time.time();answer=None
        def walk(ids,domains):
            nonlocal calls,answer
            lev=len(ids)
            for i in domains[0]:
                calls+=1
                if calls>args.cap:return False
                if calls%10000==0:print("attempts",ci,calls,"level",lev,"seconds",time.time()-start,flush=True)
                if lev+1==len(K):
                    answer=ids+[int(i)]
                    return True
                new=[]
                for r,domain in enumerate(domains[1:],start=lev+1):
                    values=ips(int(i))[domain] if lev==0 else VG[domain]@V[i]
                    nxt=domain[values==Kp[lev,r]]
                    if not len(nxt):break
                    new.append(nxt)
                if len(new)==len(domains)-1 and walk(ids+[int(i)],new):return True
            return False
        walk([],[np.arange(len(V)) for _ in range(len(K))])
        print("embedding",ci,"calls",calls,"seconds",time.time()-start,"found",answer is not None,flush=True)
        if answer is not None:
            Bp=V[answer]
            B=np.empty_like(Bp);B[order]=Bp
            assert np.array_equal(B@G@B.T,K)
            snf=smith_normal_form(sp.Matrix(B.tolist()),domain=sp.ZZ)
            invariants=[abs(int(snf[i,i]))for i in range(len(K))]
            print("invariants",invariants,flush=True)
            if invariants!=[1]*len(K):continue
            np.savez_compressed(args.output,G=G,B=B,K=K,S=np.array(row["S"],dtype=np.int64))
            print("saved",args.output,flush=True)
            return


if __name__=="__main__":main()
