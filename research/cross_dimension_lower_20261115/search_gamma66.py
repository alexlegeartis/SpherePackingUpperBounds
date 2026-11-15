"""Capped embedding search with diverse roots and dynamic smallest domains.

Unlike the earlier global attempt cap, the budget is spread over many first
vectors. Empty or capped searches make no non-embedding claim. Any successful
embedding is checked with exact matrix products and Smith invariant factors.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
from pathlib import Path
import argparse
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
    ap.add_argument("--roots",type=int,default=80)
    ap.add_argument("--per-root",type=int,default=6000)
    ap.add_argument("--seed",type=int,default=661117)
    args=ap.parse_args()
    d=np.load(args.pool)
    G=d["G"].astype(np.int64)
    W=d["V"].astype(np.int32)
    V=np.vstack((W,-W))
    VG=(V.astype(np.int64)@G)
    assert abs(VG).max()<2**20
    VG=VG.astype(np.int32)
    assert 72*int(abs(VG).max())*int(abs(V).max())<2**31
    target=json.loads((ROOT/"research/lower_bound_return_20261115/gamma72_66_abstract.json").read_text())
    K=np.array(target["normal_gram"],dtype=np.int64)
    S=np.array(target["dual_labels"],dtype=np.int64)
    rng=np.random.default_rng(args.seed)
    roots=rng.choice(len(W),min(args.roots,len(W)),replace=False)
    # Add basis rays to cover roots that random orbit sampling can miss.
    for v in np.eye(72,dtype=np.int32):
        hits=np.flatnonzero(np.all(W==v,axis=1))
        if len(hits):roots=np.r_[roots,hits[0]]
    history=[]
    found=None
    started=time.monotonic()
    for root in roots:
        cache={}
        nodes=0
        capped=False
        answer=None
        def ips(i):
            if i not in cache:
                if len(cache)>32:cache.pop(next(iter(cache)))
                values=VG@V[i]
                assert np.min(values)>=-8 and np.max(values)<=8
                cache[i]=values.astype(np.int8)
            return cache[i]
        root_ips=ips(int(root))
        initial={j:np.flatnonzero(root_ips==K[0,j]) for j in range(1,6)}
        maximum_depth=1
        def walk(assignment,domains):
            nonlocal nodes,capped,answer,maximum_depth
            maximum_depth=max(maximum_depth,len(assignment))
            if not domains:
                answer=[assignment[i] for i in range(6)]
                return True
            var=min(domains,key=lambda j:len(domains[j]))
            choices=domains[var]
            if len(choices)>1:choices=choices[rng.permutation(len(choices))]
            for v in choices:
                nodes+=1
                if nodes>args.per_root:
                    capped=True
                    return False
                new={}
                for j,domain in domains.items():
                    if j==var:continue
                    values=VG[domain]@V[v]
                    keep=domain[values==K[var,j]]
                    if not len(keep):break
                    new[j]=keep
                if len(new)==len(domains)-1 and walk({**assignment,var:int(v)},new):
                    return True
                if capped:return False
            return False
        walk({0:int(root)},initial)
        row={"root":int(root),"attempts":nodes,"capped":capped,"maximum_depth":maximum_depth,
             "found":answer is not None,"seconds_cumulative":time.monotonic()-started}
        history.append(row)
        print(json.dumps(row),flush=True)
        if answer is not None:
            B=V[answer].astype(np.int64)
            assert np.array_equal(B@G@B.T,K)
            snf=smith_normal_form(sp.Matrix(B.tolist()),domain=sp.ZZ)
            invariants=[abs(int(snf[i,i])) for i in range(6)]
            row["smith_invariants"]=invariants
            if invariants==[1]*6:
                np.savez_compressed(args.output,G=G,B=B,K=K,S=S)
                found=str(args.output)
                break
    report={"seed":args.seed,"pool_unsigned":len(W),"roots_tried":len(history),
            "per_root_cap":args.per_root,"history":history,"certificate":found,
            "nonembedding_established":False,"scope":"Finite randomized embedding discovery; successful integer witnesses require exact audit."}
    args.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+"\n")


if __name__=="__main__":main()
