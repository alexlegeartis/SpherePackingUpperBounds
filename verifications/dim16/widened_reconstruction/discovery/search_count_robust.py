"""Bounded PSD-preserving local-harmonic block reweighting; discovery only."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import linprog,minimize,minimize_scalar
from screen import Witness,triple_grid
D=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cap',type=float,default=.501)
    ap.add_argument('--rounds',type=int,default=20);ap.add_argument('--max-weight',type=float,default=30)
    ap.add_argument('--resume',action='store_true');ap.add_argument('--margin',type=float,default=1e-6)
    args=ap.parse_args();c=args.cap;w=Witness();start=time.monotonic();degree=24; margin=args.margin
    origin=np.r_[(w.L[17]@np.ones(17))**2,17*w.margins[17]]; m=sum(w.L[17+k].shape[0]+1 for k in range(17))
    Bpair=2*w.B[0,1]+w.B[1,1];Bcycle=3*w.B[1,1]
    scales=np.r_[np.maximum(origin,1e-8),np.ones(m-len(origin))]
    def pm(x):return w.gegen(np.asarray(x),16,degree)[:,1:]
    def parts(u,v,t):
        u,v,t=np.broadcast_arrays(u,v,t);g=t-u*v;z=(1-u*u)*(1-v*v)
        qs=[np.ones_like(u),g]
        for k in range(1,16):qs.append(((2*k+13)*g*qs[-1]-k*z*qs[-2])/(k+13))
        out=[]
        for k in range(17):
            pu=w.gegen(u,16+2*k,16-k);pv=w.gegen(v,16+2*k,16-k);L=w.L[17+k]
            out.extend([qs[k]*(pu@L.T)[:,i]*(pv@L.T)[:,i] for i in range(L.shape[0])]);out.append(qs[k]*w.margins[17+k]*np.sum(pu*pv,axis=-1))
        return np.stack(out,axis=-1)
    def cp(x):
        u,v,t=x.T
        return parts(u,v,t)+parts(u,t,v)+parts(v,t,u)
    def pp(x):
        ans=parts(x,x,np.ones_like(x));p=w.gegen(x,16,16)
        cross=(p@w.L[17].T)*(w.L[17]@np.ones(17))[None,:]
        ans[:,:17]+=2*cross
        ans[:,17]+=2*w.margins[17]*p.sum(axis=-1)
        return ans
    def cm(x):return sum(pm(x[:,i]) for i in range(3))
    pair_x=np.unique(np.r_[np.linspace(-1,c,16001),np.linspace(.49,c,5001)])
    triples=np.c_[triple_grid(c)]
    triples=np.vstack([triples,[[-1,-c,c],[c,c,c],[-.95,-.6,c]]])
    checkpoint=D/f'reweight_rows_cap{c:.6f}.npz'
    if args.resume and checkpoint.exists():
        old=np.load(checkpoint);pair_x=old['pair'];triples=old['triples']
    if args.resume:
        rng=np.random.default_rng(161130);u=rng.uniform(-1,c,10000);v=rng.uniform(-1,c,10000)
        u,v=np.minimum(u,v),np.maximum(u,v);radius=np.sqrt((1-u*u)*(1-v*v))
        lo=np.maximum(v,u*v-radius);hi=np.minimum(c,u*v+radius);ok=lo<=hi
        u,v,lo,hi=u[ok],v[ok],lo[ok],hi[ok]
        t=lo+(hi-lo)*rng.random(len(u));triples=np.vstack([triples,np.c_[u,v,t]])
    # Dense feasible antipodal boundary rows prevent missing narrow positive lobes.
    extra=[]
    for u in -1+np.geomspace(1e-9,.1,55):
      for v in np.linspace(-c,c,90):
        rr=np.sqrt((1-u*u)*(1-v*v));lo=max(v,u*v-rr);hi=min(c,u*v+rr)
        if lo<=hi:
          extra.extend([[u,v,t] for t in np.linspace(lo,hi,7)])
    triples=np.vstack([triples,np.array(extra)])
    for u in np.linspace(-1,c,121):
      vs=np.linspace(u,c,121);det=1+2*u*vs*c-u*u-vs*vs-c*c
      ok=det>=0
      triples=np.vstack([triples,np.c_[np.full(sum(ok),u),vs[ok],np.full(sum(ok),c)]])
    triples=np.unique(triples,axis=0)
    T=cp(triples);R=cm(triples);P=pp(pair_x);F=pm(pair_x)
    # Variables: scaled local block weights, global f coefficients, local R coefficients.
    objective=np.r_[origin/scales[:18],np.zeros(m-18),np.ones(degree),np.zeros(degree)]
    bounds=[(0,args.max_weight*scales[i]) for i in range(m)]+[(0,None)]*(2*degree)
    options=dict(time_limit=240,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9)
    history=[]
    for rnd in range(args.rounds):
        A=np.vstack([np.c_[P/scales,F,np.ones((len(P),degree))],np.c_[T/scales,np.zeros((len(T),degree)),R]])
        rhs=np.r_[np.full(len(P),-1-Bpair-margin),np.full(len(T),-Bcycle-margin)]
        rownorm=np.maximum(abs(A).max(axis=1),1e-12);A/=rownorm[:,None];rhs/=rownorm
        fit=linprog(objective,A_ub=A,b_ub=rhs,bounds=bounds,method='highs',options=options)
        if not fit.success:raise RuntimeError(fit.message)
        theta=fit.x[:m]/scales;a=fit.x[m:m+degree];r=fit.x[m+degree:]
        diag=sum(r)
        def cv(x):
            z=np.asarray(x)[None,:];return float(cp(z)[0]@theta+cm(z)[0]@r+Bcycle)
        def pv(x):
            z=np.array([x]);return float(1+Bpair+pp(z)[0]@theta+pm(z)[0]@a+diag)
        vals=T@theta+R@r+Bcycle;new=[];worst=float(max(vals))
        def det(x):p,q,t=x;return 1+2*p*q*t-p*p-q*q-t*t
        for idx in np.r_[np.argsort(vals)[-64:],np.random.default_rng(rnd).choice(len(triples),64,replace=False)]:
            opt=minimize(lambda x:-cv(x),triples[idx],method='SLSQP',bounds=[(-1,c)]*3,
                constraints=[{'type':'ineq','fun':det}],options={'maxiter':140,'ftol':1e-12})
            if det(opt.x)>=-1e-10:
                value=cv(opt.x);worst=max(worst,value)
                if value> -.95*margin:new.append(opt.x)
        pairvalues=1+Bpair+P@theta+F@a+diag;newp=[];pairmax=float(max(pairvalues))
        ids=np.r_[np.argmax(pairvalues),np.where((pairvalues[1:-1]>=pairvalues[:-2])&(pairvalues[1:-1]>=pairvalues[2:]))[0]+1]
        for i in ids:
            opt=minimize_scalar(lambda x:-pv(x),bounds=(pair_x[max(0,i-1)],pair_x[min(len(pair_x)-1,i+1)]),method='bounded',options={'xatol':1e-14})
            pairmax=max(pairmax,-float(opt.fun))
            if -opt.fun> -.95*margin:newp.append(opt.x)
        row=dict(round=rnd,certified=False,cap=c,margin=margin,radius_ratio=1/np.sqrt(2*(1-c)),objective=1+w.B[0,0]+fit.fun,
            theta=theta.tolist(),global_coefficients=a.tolist(),local_angular_coefficients=r.tolist(),added_diagonal=float(diag),
            pair_sampled_maximum=float(max(pairvalues)),cyclic_sampled_maximum=float(max(vals)),
            pair_separated_maximum=pairmax,cyclic_separated_maximum=worst,new_pair_rows=len(newp),new_cyclic_rows=len(new),
            pair_rows=len(pair_x),cyclic_rows=len(triples),max_scaled_constraint_residual=float(max(A@fit.x-rhs)),elapsed=time.monotonic()-start)
        history.append(row)
        (D/f'reweight_cap{c:.6f}.json').write_text(json.dumps(dict(scope='Numerical finite-row block reweighting; all PSD factors unchanged, nonnegative weights. No continuum pair/cyclic certificate or source SOS replay. Not a count or density bound.',degree=degree,max_weight=args.max_weight,history=history),indent=2)+'\n')
        print({k:v for k,v in row.items() if k not in ['theta','global_coefficients','local_angular_coefficients']},flush=True)
        if not new and not newp:break
        if new:
            e=np.array(new);triples=np.vstack([triples,e]);T=np.vstack([T,cp(e)]);R=np.vstack([R,cm(e)])
        if newp:
            pair_x=np.unique(np.r_[pair_x,newp]);P=pp(pair_x);F=pm(pair_x)
        np.savez_compressed(checkpoint,pair=pair_x,triples=triples)
if __name__=='__main__':main()
