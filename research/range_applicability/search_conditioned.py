"""Uncertified discovery of solid-harmonic cyclic density certificates.

The PSD kernel is supported on {0} union [r,R].  It uses radial hats in
squared radius and solid spherical harmonics.  Global one-variable signs
are constrained by Bernstein coefficients and Taylor tails.  Cyclic and
far constraints are sampled and separated numerically; this is NOT proof.
"""
import argparse, itertools, json, math, time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.special import eval_genlaguerre, gammaln, factorial, comb
from scipy.optimize import minimize
import cvxpy as cp


def solve(n, degree=31, nodes=5, angular=16, ratio=1.2, subdivisions=2,
          rounds=5, pairwise=False, output=None):
    a=n/2; r2=n/8; R2=r2*ratio
    radial=np.linspace(r2,R2,nodes); scale=np.exp(-np.pi*radial/2)
    ell_scale=np.array([scale*(r2/radial)**(ell/2) for ell in range(angular+1)])
    norms=np.exp(gammaln(np.arange(degree+1)+a)-gammaln(a)-gammaln(np.arange(degree+1)+1))
    coeff=cp.Variable(degree+1)
    blocks=[cp.Variable((nodes+1,nodes+1),symmetric=True)]+[cp.Variable((nodes,nodes),symmetric=True) for _ in range(angular)]
    core=[blocks[0][1:,1:]]+blocks[1:]
    vector=cp.hstack([cp.reshape(B,(nodes*nodes,),order='C') for B in core])
    cons=[B>>1e-4*np.eye(B.shape[0]) for B in blocks]
    signs=(-1.)**np.arange(degree+1)
    cons += [signs@coeff==1]
    def lag(x): return np.array([eval_genlaguerre(k,a-1,2*np.asarray(x))/norms[k] for k in range(degree+1)]).T
    def jets(x):
        J=np.zeros((degree+1,degree+1))
        for j in range(degree+1):
            for k in range(j,degree+1):
                J[j,k]=(-2.)**j*eval_genlaguerre(k-j,a+j-1,2*x)/norms[k]/factorial(j)
        return J
    def bern(l,h):
        C=np.zeros((degree+1,degree+1))
        for i in range(degree+1):
            for j in range(i+1): C[i,j]=comb(i,j)/comb(degree,j)*h**j
        return C@jets(l)
    def global_matrix(lo):
        tail=4*degree+n+20
        grid=np.linspace(lo,tail,math.ceil((tail-lo)/0.5)+1)
        M=np.vstack([bern(l,h-l) for l,h in zip(grid[:-1],grid[1:])]+[jets(tail)])
        return M/np.maximum(np.max(np.abs(M),axis=1,keepdims=True),1e-200)
    cons += [global_matrix(0)@cp.multiply(signs,coeff)>=1e-10]
    cons += [global_matrix(np.pi*R2-1e-10)@coeff<=-1e-10]
    def weights(u,ell=0):
        z=(u-r2)/(R2-r2)*(nodes-1); i=min(nodes-2,max(0,int(z)))
        w=np.zeros(nodes); w[i]=1-z+i; w[i+1]=z-i
        return w*ell_scale[ell]
    def harmonic(u,v,w):
        """Q_l=(sqrt(uv)/r2)^l P_l((u+v-w)/(2sqrt(uv)))."""
        x=(u+v-w)/(2*r2); uv=u*v/r2**2
        q=[1.,x]
        for ell in range(1,angular):
            q.append(((2*ell+n-2)*x*q[-1]-ell*uv*q[-2])/(ell+n-2))
        return np.asarray(q)
    def kernel_row(u,v,w):
        q=harmonic(u,v,w)
        return np.concatenate([q[ell]*np.outer(weights(u,ell),weights(v,ell)).ravel() for ell in range(angular+1)])
    def cyclic_row(u,v,w):
        return (kernel_row(u,v,w)+kernel_row(u,w,v)+kernel_row(v,w,u))*np.exp(np.pi*(u+v+w)/3)
    def far_row(u,v,z):
        cutoff=r2 if pairwise else R2
        t=-1+z*(1+(u+v-cutoff)/(2*np.sqrt(u*v)))
        w=u+v-2*np.sqrt(u*v)*t
        return kernel_row(u,v,w)*np.exp(np.pi*(u+v)/2)
    def pair_con(u):
        w=weights(u)*np.exp(np.pi*u/2)
        diag=sum(cp.sum(cp.multiply(np.outer(weights(u,ell),weights(u,ell))*np.exp(np.pi*u)*(u/r2)**ell,B)) for ell,B in enumerate(core))
        cross=2*np.exp(np.pi*u/2)*(w@blocks[0][1:,0])
        return lag(np.pi*u)@coeff+cross+diag<=-1e-6
    cons += [pair_con(u) for u in np.linspace(r2,R2,16*(nodes-1)+1)]
    fine=np.linspace(r2,R2,subdivisions*(nodes-1)+1)
    far_points=[(u,v,z) for i,u in enumerate(fine) for v in fine[i:] for z in np.linspace(0,1,41)]
    cyc_points=list(itertools.combinations_with_replacement(fine,3)) if not pairwise else []
    far_mat=sparse.csr_matrix(np.array([far_row(*p) for p in far_points]))
    cons += [far_mat@vector<=-1e-6]
    if cyc_points: cons += [sparse.csr_matrix(np.array([cyclic_row(*p) for p in cyc_points]))@vector<=-1e-6]
    objective=cp.sum(coeff)+blocks[0][0,0]
    T=np.polynomial.Polynomial([0.,1.]); gegen=[np.polynomial.Polynomial([1.]),T]
    for ell in range(1,angular):
        gegen.append(((2*ell+n-2)*T*gegen[-1]-ell*gegen[-2])/(ell+n-2))
    start=time.time(); history=[]; rng=np.random.default_rng(20261111+n)
    for rnd in range(rounds):
        p=cp.Problem(cp.Minimize(objective),cons)
        p.solve(solver='CLARABEL',max_iter=300,tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
        if coeff.value is None: break
        val=np.asarray(vector.value).ravel(); added=0; worst={}
        for kind,func,bounds in [('far',far_row,[(r2,R2),(r2,R2),(0,1)])]+([] if pairwise else [('cyclic',cyclic_row,[(r2,R2)]*3)]):
            # All small-cell corners plus random starts. Numerical separation only.
            starts=np.column_stack([rng.uniform(lo,hi,500) for lo,hi in bounds])
            if kind=='cyclic': starts=np.vstack([starts,np.asarray(cyc_points)])
            else:
                starts=np.vstack([starts,np.asarray(far_points)[::4]])
                # Isolate every real angular critical point at every radial
                # grid pair. Top-k random/grid starts alone miss endpoint
                # radial fibers whose worst angular lobe is very narrow.
                extrema=[]
                Bvals=[np.asarray(B.value) for B in core]
                for ii,u in enumerate(fine):
                    for v in fine[ii:]:
                        outer=np.outer(weights(u),weights(v))
                        cs=[np.sum(np.outer(weights(u,ell),weights(v,ell))*B)*(np.sqrt(u*v)/r2)**ell for ell,B in enumerate(Bvals)]
                        pol=sum((c*g for c,g in zip(cs,gegen)),np.polynomial.Polynomial([0.]))
                        cu=(u+v-(r2 if pairwise else R2))/(2*np.sqrt(u*v))
                        ts=[-1.,cu]+[float(z.real) for z in pol.deriv().roots() if abs(z.imag)<1e-7 and -1<z.real<cu]
                        extrema.extend((u,v,(t+1)/(cu+1)) for t in ts)
                starts=np.vstack([starts,np.asarray(extrema)])
                for q in extrema:
                    value=func(*q)@val
                    if value>-2e-7:
                        cons += [func(*q)@vector<=-2e-6];added+=1
            values=np.asarray([func(*q)@val for q in starts])
            candidates=starts[np.argsort(values)[-32:]]
            peaks=[]
            for q in candidates:
                opt=minimize(lambda z:-func(*z)@val,q,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-14,'gtol':1e-10,'maxiter':100})
                peaks.append((float(-opt.fun),opt.x))
            worst[kind]=max(x[0] for x in peaks)
            for value,point in peaks:
                if value>-2e-7:
                    cons += [func(*point)@vector<=-2e-6]; added+=1
        h={'round':rnd,'objective':float(p.value),'density_candidate':float(p.value)*np.pi**a/math.gamma(a+1)*(r2/4)**a,'status':str(p.status),'added':added,'worst_sampled':worst,'seconds':time.time()-start}
        history.append(h); print(h,flush=True)
        result={'dimension':n,'certified':False,'warning':'Uncertified cyclic solid-harmonic candidate with per-harmonic diagonal balancing; complete exact verification is mandatory.','pairwise_control':pairwise,'feature_psd_floor':1e-4,'minimum_numerical_feature_eigenvalue':min(float(np.linalg.eigvalsh(B.value)[0]) for B in blocks),'radial_degree':degree,'radial_nodes':nodes,'angular_degree':angular,'exclusion_squared':str(r2),'truncation_squared':str(R2),'solid_harmonic_scale_squared':str(r2),'status':str(p.status),'history':history}
        sc=np.r_[1,scale]
        result.update(coefficients=coeff.value.tolist(),matrices=[(blocks[0].value*np.outer(sc,sc)).tolist()]+[(B.value*np.outer(ell_scale[ell],ell_scale[ell])).tolist() for ell,B in enumerate(blocks[1:],start=1)],objective=float(p.value),density_candidate=h['density_candidate'],max_constraint_violation=max(float(np.max(c.violation())) for c in cons if c.value() is not None))
        if output:
            target=Path(output);temporary=target.with_suffix(target.suffix+'.tmp')
            temporary.write_text(json.dumps([result],indent=2));temporary.replace(target)
        if not added: break
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--dimension',type=int,default=17);ap.add_argument('--degree',type=int,default=31);ap.add_argument('--nodes',type=int,default=5);ap.add_argument('--angular',type=int,default=16);ap.add_argument('--ratio',type=float,default=1.2);ap.add_argument('--subdivisions',type=int,default=2);ap.add_argument('--rounds',type=int,default=5);ap.add_argument('--pairwise',action='store_true');ap.add_argument('--output',default='candidate.json')
    args=vars(ap.parse_args());args['n']=args.pop('dimension');solve(**args)
