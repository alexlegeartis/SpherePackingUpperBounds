#!/usr/bin/env python3
"""Independent integer checks of the support obstruction and search outputs.

This script does not turn unsuccessful heuristic searches into impossibility
proofs. Only the stated rank-two, length-seventeen instance is excluded exactly.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations,product
from pathlib import Path
import json
import numpy as np

P=Path(__file__).resolve().parent
MUL=np.array([[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]],dtype=np.uint8)

def support_file(name):
    data=list(map(int,(P/(name+'.supports')).read_text().split()))
    n,m=data[:2];supports=data[2:]
    assert len(supports)==len(set(supports))==m
    assert all(0<s<(1<<n) and s.bit_count() in (6,7) for s in supports)
    return n,supports

def signed_zero(cols):
    r=len(cols[0]);sums={(0,)*r}
    for col in cols:
        sums={tuple((a+sign*b)%3 for a,b in zip(x,col)) for x in sums for sign in (-1,1)}
    return (0,)*r in sums

def main():
    n,supports=support_file('n17k9')
    G=np.load(P/'n17k9_G.npy')
    assert G.shape==(9,17)
    C=np.zeros((1,17),dtype=np.uint8)
    for row in G:C=np.concatenate([C^MUL[a,row] for a in range(4)])
    weights=np.count_nonzero(C,axis=1)
    assert np.count_nonzero(weights==0)==1 and min(weights[weights>0])==7
    exact=set(np.sum((C[weights==7]!=0).astype(np.uint32)*(1<<np.arange(n,dtype=np.uint32)),axis=1).tolist())
    assert exact==set(supports) and len(exact)==408
    pair_counts=Counter()
    for s in supports:
        for pair in combinations([i for i in range(n) if s>>i&1],2):pair_counts[pair]+=1
    assert len(pair_counts)==136 and set(pair_counts.values())=={63}
    colors=[(1,0),(0,1),(1,1),(1,2)]
    allowed=[]
    for a in range(8):
        for b in range(8-a):
            for c in range(8-a-b):
                d=7-a-b-c
                counts=(a,b,c,d)
                cols=[v for count,v in zip(counts,colors) for _ in range(count)]
                if not signed_zero(cols):allowed.append(counts)
    assert {tuple(sorted(x)) for x in allowed}=={(0,0,1,6)}
    required_pairs=Fraction(408*15,63)
    assert required_pairs.denominator!=1
    result={'rank2_length17':{'status':'exact impossibility for this outer code and all nonzero ternary rank<=2 parity columns','codewords':len(C),'weight7_words':int(np.sum(weights==7)),'supports':408,'design':[17,7,63],'allowed_sorted_color_counts':[0,0,1,6],'required_total_monochromatic_coordinate_pairs':str(required_pairs)}}
    for name in ['n18k11','n19k12','n20k13']:
        f=P/(name+'_parity.json')
        if not f.exists():continue
        report=json.loads(f.read_text());n,supports=support_file(name)
        H=report['columns'];assert len(H)==n and all(any(x) for x in H)
        bad=[s for s in supports if signed_zero([H[i] for i in range(n) if s>>i&1])]
        assert len(bad)==report['uncovered_supports']
        result[name]={'status':'search did not produce a packing','exact_remaining_bad_supports':len(bad),'total_bad_supports':len(supports),'nonexistence_proved':False}
    (P/'verified_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
