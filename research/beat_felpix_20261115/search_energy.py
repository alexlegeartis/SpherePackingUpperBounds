"""Discovery only: sample the pinned D4 kernel and optimize PSD energies."""
from pathlib import Path
import sys, json, time
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import linprog

PACKAGE=Path(__file__).resolve().parents[2]/'verifications/dim4/refinement_energy20261115'
SRC=PACKAGE/'vendor/D4/verify/stage5_kernelbox/src'
sys.path.insert(0,str(SRC))
import kformula as kf
from functools import lru_cache

# Discovery acceleration only; proof code is never patched. Fixed radii and
# radial features recur at every angular sample.
kf.A=lru_cache(maxsize=None)(kf.A)
_q0=kf.q0
_qc={}
def cached_q0(r, ra):
    key=(str(r),ra)
    if key not in _qc: _qc[key]=_q0(r,ra)
    return _qc[key]
kf.q0=cached_q0
_rw=kf.radial_weight
_rc={}
def cached_rw(row,r,region,hi):
    key=(id(row),str(r),region,hi)
    if key not in _rc: _rc[key]=_rw(row,r,region,hi)
    return _rc[key]
kf.radial_weight=cached_rw
_quad=kf.quad_form
_mc={}
def cached_quad(m,fa,fb):
    key=(id(m),tuple((k,str(v)) for k,v in fa.items()),tuple((k,str(v)) for k,v in fb.items()))
    if key not in _mc: _mc[key]=_quad(m,fa,fb)
    return _mc[key]
kf.quad_form=cached_quad

BASE=Path(__file__).resolve().parent/'scratch'
BASE.mkdir(exist_ok=True)
A=Q(133,250)

def zonal(x, degree, even=True):
    x=np.asarray(x)
    u=[np.ones_like(x), 2*x]
    for d in range(2,degree+1):
        u.append(2*x*u[-1]-u[-2])
    return np.stack([u[d]/(d+1) for d in range(0,degree+1,2 if even else 1)],axis=-1)

def samples(j):
    target=BASE/f'samples_j{j}.json'
    if target.exists(): return json.loads(target.read_text())
    ts=sorted(set([Q(-1)+(1-A)*Q(k,400) for k in range(401)] +
                  [Q(-87,100)+Q(k,10000) for k in range(121)]))
    # Sample endpoints and cell midpoints. A grid can miss one-sided jumps;
    # only the full interval tree can certify the kernel domain.
    rad=sorted(set(kf.BREAKS[:j+1]+[(kf.BREAKS[k]+kf.BREAKS[k+1])/2 for k in range(j)]))
    rows=[]
    start=time.time()
    for r in rad:
        for s in rad:
            if s<r: continue
            for t in ts:
                v=kf.minusK_point(r,s,t)
                rows.append([str(r),str(s),str(t),float(v.mid())])
            print('samples',j,str(r),str(s),len(rows),round(time.time()-start,1),flush=True)
    target.write_text(json.dumps(rows))
    return rows

def solve(j, degree, rows, even=True, credits=False):
    rows=[x for x in rows if np.isfinite(x[3])]
    tfar=np.array([float(Q(x[2])) for x in rows])
    budget=np.array([x[3]*1e6 for x in rows])
    cap=1-Q(1)/(2*kf.BREAKS[j]**2)
    ti=np.linspace(0,float(A),3001) if even else np.linspace(-float(A),float(cap),5001)
    Af=zonal(tfar,degree,even)
    Ai=zonal(ti,degree,even)
    n=24 if j==1 else 26
    size=degree//2+1 if even else degree+1
    obj=np.full(size,float(n)); obj[0]-=n*n
    bounds=[(None,None)]+[(0,None)]*(size-1)
    if credits:
        credit={r['degree']:float(Q(r['beta']))*1e6 for r in json.loads((PACKAGE/'gram_credits.json').read_text()).get(str(j),[])}
        degrees=list(range(0,degree+1,2 if even else 1))
        bounds=[(None,None)]+[(-credit.get(d,0),None) for d in degrees[1:]]
    res=linprog(obj,A_ub=np.concatenate([Af,Ai]),b_ub=np.concatenate([budget-1e-5,np.full(len(ti),-1e-5)]),
                bounds=bounds,method='highs')
    if not res.success:
        print('FAILED',j,degree,res.message,flush=True); return
    out=dict(j=j,degree=degree,even=even,credits=credits,coefficients=res.x.tolist(),energy=-res.fun,
             p1=sum(res.x),c0=res.x[0],min_margin=float(min(budget-Af@res.x)),
             angular_max=float(max(Ai@res.x)),nit=res.nit)
    print(json.dumps(out),flush=True)
    (BASE/f'energy_j{j}_d{degree}_{"even" if even else "full"}{"_credits" if credits else ""}.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':
    j=int(sys.argv[1])
    data=samples(j)
    even='--full' not in sys.argv
    credits='--credits' in sys.argv
    for d in (4,6,8,10,12,14,16,18,20,24,28,32): solve(j,d,data,even,credits)
