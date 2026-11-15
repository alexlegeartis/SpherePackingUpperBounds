"""Sampled fixed-cardinality spherical SDP for far-pair energy (discovery only).

The anchor Gram representation follows the pinned Felpix spherical proof,
including its fixed-count restriction. No solver output is a certificate.
All quantities in this discovery model are scaled by 10^6.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,time
import numpy as np
import cvxpy as cp
from scipy.optimize import minimize,minimize_scalar
from scipy.sparse import csr_matrix

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

class Model:
    def __init__(self,n,s,d,L,D):
        self.n,self.s,self.d,self.L,self.D=n,s,d,L,D
        self.nr=d+1
        self.size=self.nr*self.nr
        self.free=(L+1)*self.size
        self.pstart=self.free+d+2
        self.total=self.pstart+D+1
        self.obj=np.zeros(self.total)
        self.obj[0]=1
        self.obj[self.free]=-2*(n-1)
        self.obj[self.free+1]=-(n-1)**2

    def basis(self,a):
        x=2*(np.asarray(a)+1)/(1+self.s)-1
        return np.polynomial.chebyshev.chebvander(x,self.d)

    def zonal(self,t):
        t=np.asarray(t)
        u=[np.ones_like(t),2*t]
        for k in range(2,self.D+1):u.append(2*t*u[-1]-u[-2])
        return np.stack([x/(k+1) for k,x in enumerate(u[:self.D+1])],axis=-1)

    def kernel(self,a,b,t):
        a,b,t=np.broadcast_arrays(a,b,t)
        ba,bb=self.basis(a),self.basis(b)
        outer=(ba[..., :,None]*bb[...,None,:]).reshape((-1,self.size))
        out=np.zeros((a.size,self.total))
        # The unrestricted degree-zero matrix is reconstructed from its
        # positive restriction to [1,sum T_1,...,sum T_d].
        degree0=np.zeros_like(outer).reshape((-1,self.nr,self.nr))
        degree0[:,1:,1:]=(ba[...,1:,None]*bb[...,None,1:]).reshape((-1,self.d,self.d))
        out[:,:self.size]=degree0.reshape((-1,self.size))
        out[:,self.free+1]=1
        out[:,self.free+2:self.pstart]=ba.reshape((-1,self.nr))[:,1:]+bb.reshape((-1,self.nr))[:,1:]
        x=(t-a*b).ravel();w=((1-a*a)*(1-b*b)).ravel()
        h0=np.ones(a.size);h1=x
        for ell in range(1,self.L+1):
            out[:,ell*self.size:(ell+1)*self.size]=h1[:,None]*outer
            h0,h1=h1,((2*ell+1)*x*h1-ell*w*h0)/(ell+1)
        return out

    def pair(self,a):
        a=np.asarray(a)
        out=self.kernel(a,a,np.ones_like(a))
        ba=self.basis(a)
        for i in range(1,self.nr):out[:,i]+=2*ba[:,i]
        out[:,self.free]+=2
        out[:,self.free+2:self.pstart]-=2*(self.n-1)*ba[:,1:]
        out[:,self.pstart:]-=self.zonal(a)
        return out

    def triangle(self,triples):
        a,b,t=np.asarray(triples).T
        return self.kernel(a,b,t)+self.kernel(a,t,b)+self.kernel(b,t,a)

    def p_rows(self,t):
        out=np.zeros((len(t),self.total));out[:,self.pstart:]=self.zonal(t)
        return out


def random_triples(s,count,rng):
    a=rng.uniform(-1,s,count*3);b=rng.uniform(-1,s,count*3)
    radius=np.sqrt(np.maximum(0,(1-a*a)*(1-b*b)))
    lower=np.maximum(-1,a*b-radius);upper=np.minimum(s,a*b+radius)
    okay=lower<=upper
    a,b,lower,upper=[x[okay] for x in [a,b,lower,upper]]
    t=lower+rng.uniform(0,1,len(a))*(upper-lower)
    # Include Gram-determinant boundary faces as well.
    t[::4]=lower[::4];t[1::4]=upper[1::4]
    return np.sort(np.stack([a,b,t],axis=-1),axis=1)[:count]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--N',type=int,default=24)
    ap.add_argument('--j',type=int,choices=[1],default=1)
    ap.add_argument('--degree',type=int,default=4)
    ap.add_argument('--angular-degree',type=int,default=6)
    ap.add_argument('--p-degree',type=int,default=16)
    ap.add_argument('--rounds',type=int,default=12)
    ap.add_argument('--joint-spherical',action='store_true')
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--samples',type=Path,default=HERE/'dim4_research/kernel_samples_j1.json')
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    R=Q(501473,500000) if args.j==1 else Q(102597,100000)
    s=float(1-1/(2*R*R));a=.532
    model=Model(args.N,s,args.degree,args.angular_degree,args.p_degree)
    offset=np.zeros(model.D+1);constant=0.
    if args.joint_spherical:
        origin=json.loads((ROOT/'verifications/dim4/refinement_energy20261115/vendor/D4/data/v7/origin.json').read_text())
        for row in origin['appended_spherical_local_energy_rays']:
            weight=float(Q(row['weight']))*1e6
            constant+=weight*float(Q(row['B']))
            for k,c in enumerate(row['scalar_coefficients']):offset[k]-=weight*float(Q(c))
    def offset_value(t):return model.zonal(np.asarray(t))@offset
    rng=np.random.default_rng(20261115)
    triples=random_triples(s,2000,rng)
    special=[-1,-.875,-.75,-.5,-.25,0,.25,.5,s]
    import itertools
    boundary=[z for z in itertools.combinations_with_replacement(special,3)
              if 1+2*np.prod(z)-sum(x*x for x in z)>=-1e-14]
    triples=np.concatenate([triples,np.array(boundary)])
    pair=np.linspace(-1,s,1501)
    central=np.linspace(-a,s,2001)
    data=json.loads(args.samples.read_text())
    data=[x for x in data if np.isfinite(x[3])]
    far=np.array([float(Q(x[2])) for x in data])
    budget=np.array([x[3]*1e6 for x in data])
    matrices=[cp.Variable((model.nr,model.nr),PSD=True) for _ in range(model.L+1)]
    free=cp.Variable(model.d+2)
    p=cp.Variable(model.D+1)
    variables=cp.hstack([cp.vec(m,order='C') for m in matrices]+[free,p])
    start=time.time()
    for round_index in range(args.rounds):
        rows=np.vstack([model.triangle(triples),model.pair(pair),model.p_rows(central),model.p_rows(far)])
        rhs=np.r_[np.full(len(triples),-1e-4),offset_value(pair)-1e-4,
                  np.full(len(central),-1e-4),budget-1e-4]
        scale=np.maximum(np.max(abs(rows),axis=1),1)
        constraints=[csr_matrix(rows/scale[:,None])@variables<=rhs/scale]
        constraints += [cp.trace(m)<=1e5 for m in matrices]
        constraints += [cp.abs(free)<=1e4,cp.abs(p)<=1e4]
        problem=cp.Problem(cp.Minimize(model.obj@variables),constraints)
        attempts=[]
        accepted=False
        for tolerance in [1e-6,1e-5]:
            try:
                problem.solve(solver='CLARABEL',max_iter=250,max_threads=1,
                              tol_gap_abs=tolerance,tol_feas=tolerance,tol_gap_rel=tolerance)
                attempts.append(dict(solver='CLARABEL',tolerance=tolerance,status=problem.status))
                if problem.status in [cp.OPTIMAL,cp.OPTIMAL_INACCURATE] and variables.value is not None:
                    accepted=True
                    break
            except cp.error.SolverError as exc:attempts.append(dict(solver='CLARABEL',error=str(exc)))
        if not accepted:
            problem.solve(solver='SCS',eps=1e-5,max_iters=50000)
            attempts.append(dict(solver='SCS',status=problem.status))
        if variables.value is None:raise RuntimeError('no solution: '+problem.status)
        values=np.asarray(variables.value).ravel()
        report=dict(status=problem.status,certified=False,N=args.N,j=args.j,s=s,radial_degree=model.d,
                    angular_degree=model.L,p_degree=model.D,round=round_index,
                    kappa_micro=float(model.obj@values),energy_micro=float(args.N*(constant+sum(offset)-(model.obj@values))),
                    joint_spherical=args.joint_spherical,
                    max_sample_residual=float(np.max(rows@values-rhs)),seconds=time.time()-start,solver_attempts=attempts)
        np.savez(args.out/'proposal.npz',matrices=np.stack([m.value for m in matrices]),free=free.value,
                 p=p.value,values=values,N=args.N,j=args.j,s=s,radial_degree=model.d,
                 angular_degree=model.L,p_degree=model.D)
        (args.out/'candidate.json').write_text(json.dumps(report,indent=2)+'\n')
        # Random and numerical separation is discovery; exact Bernstein and
        # kernel-box replays are required before any density claim.
        tests=random_triples(s,20000,rng)
        residual=model.triangle(tests)@values
        new_tri=[]
        for index in np.argsort(residual)[-20:]:
            z=tests[index]
            def objective(x):return -float(model.triangle(np.array([x]))@values)
            opt=minimize(objective,z,method='SLSQP',bounds=[(-1,s)]*3,
                         constraints=[dict(type='ineq',fun=lambda x:1+2*np.prod(x)-np.dot(x,x))],
                         options={'maxiter':150,'ftol':1e-8})
            if opt.success and -opt.fun>1e-5 and 1+2*np.prod(opt.x)-np.dot(opt.x,opt.x)>=-1e-8:
                new_tri.append(np.sort(opt.x))
        pairgrid=np.linspace(-1,s,6001)
        pr=model.pair(pairgrid)@values-offset_value(pairgrid)
        new_pair=[]
        for index in range(1,len(pr)-1):
            if pr[index]>pr[index-1] and pr[index]>pr[index+1]:
                opt=minimize_scalar(lambda x:-float(model.pair(np.array([x]))@values-offset_value(np.array([x]))),
                                    bounds=(pairgrid[index-1],pairgrid[index+1]),method='bounded')
                if -opt.fun>1e-5:new_pair.append(opt.x)
        report.update(new_triangle_cuts=len(new_tri),new_pair_cuts=len(new_pair),
                      max_random_triangle=float(max(residual)),max_dense_pair=float(max(pr)))
        print(json.dumps(report),flush=True)
        (args.out/f'round{round_index:02}.json').write_text(json.dumps(report,indent=2)+'\n')
        if not new_tri and not new_pair:break
        if new_tri:triples=np.concatenate([triples,np.array(new_tri)])
        if new_pair:pair=np.r_[pair,new_pair]


if __name__=='__main__':main()
