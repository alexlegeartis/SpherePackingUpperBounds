"""Enrich the minimum shell using the six published Gamma72 automorphisms.

Every generator identity is exact. The resulting finite pool is for discovery,
not a claim to enumerate the full 6,218,175,600-vector minimum shell.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
from pathlib import Path
import argparse
import json
import re
import time
import numpy as np
from scipy.sparse import csr_matrix

ROOT=Path(__file__).resolve().parents[2]
PARENT=ROOT/"research/novel_routes_20261113/higher_parent"


def generators(G):
    text=(PARENT/"neb72.html").read_text()
    body=text.split('<a NAME="SUBGROUP_GENERATORS"',1)[1].split('</a><br>',1)[1].split('<a NAME=',1)[0]
    nums=list(map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',body))))
    p=1
    result=[]
    for _ in range(nums[0]):
        n=nums[p];p+=1
        assert n==72
        A=np.array(nums[p:p+n*n],dtype=np.int64).reshape(n,n);p+=n*n
        assert np.array_equal(A@G@A.T,G)
        result.append(csr_matrix(A.T))
    assert p==len(nums)
    return result


def canonical(V):
    V=V.copy()
    V*=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)])[:,None]
    return np.unique(V,axis=0)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--size",type=int,default=450000)
    ap.add_argument("--rounds",type=int,default=3)
    ap.add_argument("--seed",type=int,default=661116)
    args=ap.parse_args()
    d=np.load(args.input);G=d["G"].astype(np.int64);V=d["V"].astype(np.int16)
    assert np.array_equal(G,np.load(PARENT/"G72.npy"))
    assert np.all(np.einsum('vi,vi->v',V.astype(np.int64)@G,V.astype(np.int64))==8)
    gens=generators(G)
    rng=np.random.default_rng(args.seed)
    start=time.monotonic()
    history=[]
    for step in range(args.rounds):
        selected=V[rng.choice(len(V),min(100000,len(V)),replace=False)].astype(np.int64)
        pieces=[V]
        for A in gens:
            W=np.asarray((A@selected.T).T)
            assert abs(W).max()<32767
            pieces.append(W.astype(np.int16))
        V=canonical(np.vstack(pieces))
        if len(V)>args.size:
            V=V[rng.choice(len(V),args.size,replace=False)]
        row={"round":step,"pool_unsigned":len(V),"seconds":time.monotonic()-start}
        history.append(row);print(json.dumps(row),flush=True)
    assert np.all(np.einsum('vi,vi->v',V.astype(np.int64)@G,V.astype(np.int64))==8)
    np.savez_compressed(args.output,G=G,V=V)
    args.output.with_suffix('.json').write_text(json.dumps({"seed":args.seed,"history":history,
        "source_automorphism_identities_exact":True,"all_returned_norms_exact":8,
        "full_shell_enumerated":False},indent=2)+"\n")


if __name__=="__main__":main()
