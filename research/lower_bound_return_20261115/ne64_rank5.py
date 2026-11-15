"""Reproduce discovery of the six-point 59D Ne64 antipode.

Floating scores and a finite pool select candidates only. The final exact proof
is in constructions/dim59; this program does not establish global optimality.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
from pathlib import Path
import argparse
import itertools
import json
import math
import numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/"novel_routes_20261113/dual64"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--cap",type=int,default=2000000)
    args=ap.parse_args()
    pairs=list(itertools.combinations(range(5),2))
    W=np.array(list(itertools.product(range(-1,2),repeat=5)),dtype=np.int64)
    W=W[np.any(W,axis=1)]
    candidates=[]
    for weights in itertools.product((-3,-2),repeat=10):
        K=12*np.eye(5,dtype=np.int64)
        for (i,j),v in zip(pairs,weights):K[i,j]=K[j,i]=v
        if np.einsum('vi,ij,vj->v',W,K,W).min()<12:continue
        if np.linalg.eigvalsh(K)[0]<=0:continue
        H=3*np.linalg.inv(K)
        beta=max(H.diagonal().max(),(H.diagonal()[:,None]+H.diagonal()[None,:]-2*H).max())
        score=6*((12-beta)/4)**29.5/(3**16*math.sqrt(np.linalg.det(K)/3**5))
        candidates.append((score,K))
    candidates.sort(key=lambda item:-item[0])
    a=np.load(SOURCE/"reduced_dual.npz")
    V=np.load(SOURCE/"short_pool.npz")["V"].astype(np.int64)
    V=np.vstack((V,-V));VR=V@a["R"]
    K=candidates[0][1]
    assert np.all(np.einsum('vi,ij,vj->v',V,a["R"],V)==12)
    attempts=0;answer=None
    def walk(ids,domains):
        nonlocal attempts,answer
        lev=len(ids)
        for i in domains[0]:
            attempts+=1
            if attempts>args.cap:return False
            if lev+1==5:answer=ids+[int(i)];return True
            new=[]
            for r,domain in enumerate(domains[1:],start=lev+1):
                nxt=domain[(VR[domain]@V[i])==K[lev,r]]
                if not len(nxt):break
                new.append(nxt)
            if len(new)==len(domains)-1 and walk(ids+[int(i)],new):return True
        return False
    walk([],[np.arange(len(V))for _ in range(5)])
    if answer is None:raise SystemExit("No embedding within the finite discovery cap; no proof claim.")
    A=V[answer]@a["T"]
    assert np.array_equal(A@a["H"]@A.T,K)
    sf=smith_normal_form(Matrix(A.tolist()),domain=ZZ)
    assert [abs(int(sf[i,i]))for i in range(5)]==[1]*5
    S=np.vstack((np.zeros((1,5),dtype=np.int64),np.eye(5,dtype=np.int64)))
    np.savez_compressed(args.output,G=a["G"],H=a["H"],A=A,Knum=K,S=S)
    print(json.dumps({"embedding_found":True,"attempts":attempts,"pool_unsigned":len(V)//2,
                      "discovery_score":candidates[0][0],"normal_gram_numerator":K.tolist(),
                      "proof":"Run constructions/dim59/verify.py and verify_independent.py."},indent=2))


if __name__=="__main__":main()
