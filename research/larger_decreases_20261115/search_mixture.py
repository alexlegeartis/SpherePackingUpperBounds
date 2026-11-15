"""LP for a positive sum of profile kernels requiring individual certification.

Discovery only. Each component needs a complete angular replay and the joint
radial witness needs global Fourier, tail and pair verification.
"""
from pathlib import Path
from fractions import Fraction as Q
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,json,os,sys,subprocess,time,math
import numpy as np
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
from flint import arb,arb_poly,ctx
from common import source

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('dimension',type=int)
    ap.add_argument('--degree',type=int,default=95)
    ap.add_argument('--cut',type=float,default=60)
    ap.add_argument('--jobs',type=int,default=2)
    args=ap.parse_args();n=args.dimension;d=args.degree;s=Q(n,8)
    profiles=[(0,20),(0,40),(-32,56),(-96,120),(-192,216),(-384,408),
              (-96,104),(-96,112),(-96,128),(-96,144),(-96,192),(-32,40),(-32,64)]
    folder=HERE/f'mixture_n{n}_d{d}';folder.mkdir(exist_ok=True)
    def angular(ab):
        a,b=ab;tag=f'library_n{n}_a{a}_b{b}'
        dest=HERE/tag
        if not (dest/'angular.json').exists():
            r=subprocess.run([sys.executable,str(HERE/'search_angular.py'),str(n),f'--a={a}',f'--b={b}',
                              '--rounds','6','--tag',tag],cwd=ROOT,env=ENV,capture_output=True,text=True)
            (folder/(tag+'.log')).write_text(r.stdout+r.stderr)
            if r.returncode:
                print('ANGULAR_SEARCH_REJECTED',tag,flush=True)
                return None
        angular=json.loads((dest/'angular.json').read_text())
        return dict(a=a,b=b,B=float(Q(angular['B_exact'])),angular=str((dest/'angular.json').relative_to(ROOT)))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        library=[x for x in pool.map(angular,profiles) if x is not None]
    (folder/'library.json').write_text(json.dumps(library,indent=2)+'\n')
    ctx.prec=768
    p,h,_=source(n,d)
    AP=arb_poly([arb(x) for x in p.coeffs()]);AH=arb_poly([arb(x) for x in h.coeffs()])
    norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1))
    sign=(-1.)**np.arange(d+1)
    def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*np.asarray(x))/norm[k] for k in range(d+1)]).T
    def ev(poly,x):return np.array([float(poly(arb(float(z))).mid()) for z in x])
    limit=max(300.,4*d+n+20.)
    grids={'pair':np.unique(np.r_[np.linspace(float(s)*np.pi,args.cut,3501),float(s)*np.pi+np.geomspace(1e-8,3,501)]),
           'fourier':np.r_[np.linspace(0,limit,6001),np.geomspace(limit,2500,501)],
           'physical':np.r_[np.linspace(args.cut,limit,5001),np.geomspace(limit,2500,501)]}
    blocks=[]
    for kind,x in grids.items():
        if kind=='pair':
            v=x/(np.pi*float(s))
            hh=np.stack([np.exp(-z['a']*np.log(v)-z['b']*(v-1)) for z in library],axis=-1)
            allowance=np.exp(x)[:,None]*(2*hh-hh*hh)/np.array([z['B'] for z in library])[None,:]
            row=np.c_[ev(AP,x),lag(x),-allowance]
        elif kind=='fourier':row=np.c_[-ev(AH,x),-lag(x)*sign,np.zeros((len(x),len(library)))]
        else:row=np.c_[ev(AP,x),lag(x),np.zeros((len(x),len(library)))]
        row/=np.maximum(abs(row).max(axis=1,keepdims=True),1e-200)
        blocks.append(row)
    A=np.vstack(blocks)
    obj=np.r_[float(p[0]),np.ones(d+1+len(library))]
    eq=np.r_[float(h[0]),sign,np.zeros(len(library))]
    bounds=[(0,None)]+[(None,None)]*(d+1)+[(0,None)]*len(library)
    start=time.time()
    res=linprog(obj,A_ub=A,b_ub=np.full(len(A),-1e-7),A_eq=eq[None,:],b_eq=[1],bounds=bounds,
                method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'time_limit':120})
    if not res.success:raise RuntimeError(res.message)
    density=res.fun*np.pi**(n/2)*(float(s)/4)**(n/2)/math.gamma(n/2+1)
    old=next(z for z in json.loads((ROOT/'RESULTS.json').read_text()) if z['dimension']==n)
    out=dict(dimension=n,degree=d,s_exact=str(s),cut=args.cut,library=library,x=res.x.tolist(),
             objective=res.fun,density_candidate=density,previous_strict=old['strict_density_upper_bound'],
             numerical_decrease_percent=100*(1-density/float(old['strict_density_upper_bound'])),
             active_components=[dict(**library[i],weight=float(w)) for i,w in enumerate(res.x[d+2:]) if w>1e-5],
             max_scaled_residual=float(max(A@res.x)),seconds=time.time()-start,certified=False)
    (folder/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['library','x']}),flush=True)

if __name__=='__main__':main()
