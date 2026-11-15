"""Reoptimize the scalar Fourier function against the frozen 4D kernel.

K_new=u*K_old. P_new<=u*P_old on 2<=x<=49/8 and P_new<=0 above
49/8 preserve every pair inequality. H_new>=0 and H_new(0)=1 leave
density <= (P_new(0)+u*(A-floor))/8. Discovery only.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,math,time
import numpy as np
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--degree',type=int,default=47)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--margin',type=float,default=1e-9)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    base=ROOT/'verifications/dim4/refinement_energy20261115'
    census=json.loads((base/'vendor/D4/data/census_floors.json').read_text())
    certificate=json.loads((base/'certificate.json').read_text())
    scalar=json.loads((base/'vendor/D4/data/v7/scalar.json').read_text())
    old=list(map(Q,scalar['F']))
    assert Q(scalar['alpha'])==2
    old_upper=Q(certificate['exact_density_upper_bound'])
    # P_old(0)+A-floor=8*old_upper, with the broad-graph saving included.
    kernel_cost=8*old_upper-old[0]
    d=args.degree;sign=(-1.)**np.arange(d+1)
    def lag(x):return np.array([eval_genlaguerre(k,1,2*np.asarray(x))/(k+1) for k in range(d+1)]).T
    def original(x):return np.polynomial.polynomial.polyval(x,np.array(list(map(float,old))))
    cutoff=49/8
    grids=dict(pair=np.linspace(2,cutoff,3501),physical=np.r_[np.linspace(cutoff,4*d+50,6501),np.geomspace(4*d+50,2000,501)],
               fourier=np.r_[np.linspace(0,4*d+50,7501),np.geomspace(4*d+50,2000,501)])
    cutsfile=args.out/'cuts.json'
    if cutsfile.exists():
        for k,v in json.loads(cutsfile.read_text()).items():grids[k]=np.r_[grids[k],v]
    rows=[]
    for kind,x in grids.items():
        if kind=='pair':r=np.c_[lag(x),-original(x)]
        elif kind=='physical':r=np.c_[lag(x),np.zeros(len(x))]
        else:r=np.c_[-lag(x)*sign,np.zeros(len(x))]
        r/=np.maximum(abs(r).max(axis=1,keepdims=True),1e-200)
        rows.append(r)
    A=np.vstack(rows)
    objective=np.r_[np.ones(d+1),float(kernel_cost)]
    equality=np.r_[sign,0]
    start=time.time();attempts=[]
    for method in ['highs','highs-ipm','highs-ds']:
        result=linprog(objective,A_ub=A,b_ub=np.full(len(A),-args.margin),
                       A_eq=equality[None,:],b_eq=[1],bounds=[(None,None)]*(d+1)+[(0,None)],method=method,
                       options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'time_limit':120})
        attempts.append(dict(method=method,status=result.message))
        if result.success:break
    if not result.success:raise RuntimeError(result.message)
    out=dict(dimension=4,degree=d,coefficients=result.x[:-1].tolist(),kernel_scale=result.x[-1],
             kernel_cost_exact=str(kernel_cost),density_candidate=result.fun/8,
             previous_strict='0.629520189662',decrease_percent=100*(1-result.fun/8/float(Q('0.629520189662'))),
             max_scaled_residual=float(max(A@result.x)),seconds=time.time()-start,solver_attempts=attempts,certified=False)
    (args.out/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['coefficients','kernel_cost_exact']}),flush=True)

if __name__=='__main__':main()
