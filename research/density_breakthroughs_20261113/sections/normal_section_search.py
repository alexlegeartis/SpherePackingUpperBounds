"""Bounded discovery search for codimension-two sections avoiding norm-eight.

This is only a necessary-condition search.  A successful normal pair would still
need primitive saturation and an independent restricted minimum calculation.
"""
import json, time, sys
from pathlib import Path
import numpy as np

P=Path(__file__).resolve().parent
d=np.load(P/'normal_contact_coordinates.npz')
V=d['V'].astype(np.int16); H=d['H']; n=H.shape[0]
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 7892)
seconds=float(sys.argv[1]) if len(sys.argv)>1 else 90.
hard=len(sys.argv)>3 and sys.argv[3]=='hard'
output=P/('normal_section_search_hard.json' if hard else 'normal_section_search.json')
limit=2785.573433034726
start=time.monotonic(); best={}; records=[]; trial=0

def measure(C):
    M=C@H@C.T
    return float(np.linalg.det(M)), int(np.sum(np.all(V@C.T==0,axis=1)))

while time.monotonic()-start<seconds:
    C=rng.choice([-1,0,1],size=(2,n),p=[.2,.6,.2]).astype(np.int64)
    if hard and measure(C)[0]>limit:continue
    lam=float(rng.choice([1,3,8,15,30,60]))
    old=10**9; stagnation=0
    for sweep in range(150):
        if time.monotonic()-start>=seconds: break
        for j in [0,1]:
            c=C[j].copy(); other=C[1-j]
            q=H@other; norm=other@q
            if norm<1e-8: break
            D=norm*H-np.outer(q,q)
            A=V[(V@other)==0].astype(np.int64)
            y=A@c; Dc=D@c; det=c@Dc
            for i in rng.permutation(n):
                shifts=np.arange(-3,4,dtype=np.int64)
                zeros=(y[:,None]+A[:,i,None]*shifts[None,:]==0).sum(0)
                dets=det+2*shifts*Dc[i]+shifts**2*D[i,i]
                costs=zeros+lam*dets/limit
                costs[dets<1e-7]=np.inf
                if hard:costs[dets>limit]=np.inf
                k=int(np.argmin(costs)); step=shifts[k]
                if step:
                    c[i]+=step; y+=step*A[:,i]
                    det=dets[k]; Dc+=step*D[:,i]
            C[j]=c
        det,zeros=measure(C)
        # Keep the lowest number of surviving contacts at several volume caps.
        changed=False
        for cap in [limit,2*limit,4*limit,10*limit]:
            key=str(round(cap,5))
            prev=best.get(key)
            if 1e-7<det<=cap and (prev is None or (zeros,det)<(prev['zeros'],prev['determinant'])):
                row=dict(trial=trial,sweep=sweep,zeros=zeros,determinant=det,normal_coefficients=C.tolist(),elapsed=time.monotonic()-start,penalty=lam)
                best[key]=row; records.append(dict(cap=cap,**row));changed=True
                print(json.dumps({k:v for k,v in records[-1].items() if k!='normal_coefficients'}),flush=True)
        if changed:
            output.write_text(json.dumps(dict(threshold=limit,hard=hard,best=best,records=records),indent=2))
        if zeros==0 and det<limit:
            np.savez(P/'shell_avoiding_normals_candidate.npz',C=C,H=H)
            print('NECESSARY CONDITION SUCCESS',flush=True);sys.exit(0)
        obj=zeros+lam*det/limit
        if obj>=old-1e-6:stagnation+=1
        else:stagnation=0
        old=obj
        if stagnation>=3:
            for _ in range(2+stagnation//4):
                C[rng.integers(2),rng.integers(n)]+=rng.choice([-1,1])
        if stagnation>=15:break
    trial+=1
print(json.dumps(dict(done=True,trials=trial,elapsed=time.monotonic()-start)),flush=True)
