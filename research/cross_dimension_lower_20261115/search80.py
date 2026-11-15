"""Find small antipode normals in published extremal L80/M80 parents.

Floating scores select actual integer embeddings. Final proof is in verify80.py.
No worldwide record conclusion is drawn from the sparse published tables.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
from pathlib import Path
import argparse
import itertools
import json
import time
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ROOT=Path(__file__).resolve().parents[2]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pool",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--parent",default="L80")
    ap.add_argument("--cap",type=int,default=100000)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    data=np.load(args.pool);G=data["G"].astype(np.int64);V=data["V"].astype(np.int64)
    V=np.vstack((V,-V));VG=V@G
    assert np.all(np.einsum('vi,vi->v',VG,V)==8)
    rows=[]
    for k in range(1,5):
        n=80-k
        S=np.vstack((np.zeros((1,k),dtype=np.int64),np.eye(k,dtype=np.int64)))
        calls=0;answer=None
        if k==1:
            i=int(np.flatnonzero(G.diagonal()==10)[0]);B=np.eye(80,dtype=np.int64)[[i]]
            K=B@G@B.T
        else:
            K=np.load(ROOT/f"research/novel_routes_20261113/higher_parent/gamma72_simplex_rank{k}.npz")["K"]
            def walk(ids,domains):
                nonlocal calls,answer
                lev=len(ids)
                for i in domains[0]:
                    calls+=1
                    if calls>args.cap:return False
                    if lev+1==k:answer=ids+[int(i)];return True
                    new=[]
                    for r,domain in enumerate(domains[1:],start=lev+1):
                        keep=domain[(VG[domain]@V[i])==K[lev,r]]
                        if not len(keep):break
                        new.append(keep)
                    if len(new)==len(domains)-1 and walk(ids+[int(i)],new):return True
                return False
            walk([],[np.arange(len(V)) for _ in range(k)])
            if answer is None:
                rows.append({"dimension":n,"embedding_found":False,"attempts":calls})
                print(json.dumps(rows[-1]),flush=True)
                continue
            B=V[answer]
        assert np.array_equal(B@G@B.T,K)
        snf=smith_normal_form(sp.Matrix(B.tolist()),domain=sp.ZZ)
        assert all(abs(snf[i,i])==1 for i in range(k))
        H=sp.Matrix(K.tolist()).inv()
        beta=sp.Rational(0)
        for x,y in itertools.combinations(S,2):
            delta=sp.Matrix([(x-y).tolist()])
            beta=max(beta,(delta*H*delta.T)[0])
        d2=sp.Rational(len(S)**2,sp.Matrix(K.tolist()).det())*((8-beta)/4)**n
        name=f"{args.parent}_n{n}.npz"
        np.savez_compressed(args.output/name,G=G,B=B,K=K,S=S)
        rows.append({"parent":args.parent,"dimension":n,"certificate":name,"attempts":calls,
                     "cosets":len(S),"beta":str(beta),"density_squared":str(d2),
                     "center_density":float(sp.sqrt(d2)),"exact_audit_pending":True,
                     "published_record_improvement_established":False})
        print(json.dumps(rows[-1]),flush=True)
    (args.output/(args.parent+"_search.json")).write_text(json.dumps(rows,indent=2)+"\n")


if __name__=="__main__":main()
