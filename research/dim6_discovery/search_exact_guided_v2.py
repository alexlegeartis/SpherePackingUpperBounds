"""Uncertified discovery with Chebyshev-radius cyclic kernels and count cuts.

The PSD kernel uses ordinary normalized Gegenbauer angular factors and
Chebyshev polynomials in radius on a compact annulus, plus an isolated origin.
Selected variants freeze the Fourier polynomial, preserve separation cuts,
or add pair/antipodal separation. Numerical solver output is not a bound.
Complete exact Fourier, PSD, pair, cyclic and far verification is mandatory.
"""
import argparse, itertools, json, math, time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.special import eval_genlaguerre, gammaln, factorial, comb
from scipy.optimize import minimize, minimize_scalar
import cvxpy as cp


def solve(n, degree=31, nodes=5, angular=16, ratio=1.2, subdivisions=2,
          rounds=5, pairwise=False, output=None, exclusion=1.0, counts=False, seed=None):
    a=n/2; r2=exclusion; R2=r2*ratio
    count_table={5:[(.003012,45),(.008096,46),(.013258,47),(.017445,48)],6:[(.002005,79),(.004023,80),(.006053,81),(.008096,82)],7:[(.001000,136),(.002005,137),(.003012,138),(.004023,139),(.005036,140)],9:[(.003012,373),(.029232,457),(.030324,459),(.031420,464),(.032519,468),(.033621,473)]}
    inputs=count_table.get(n,[]) if counts else []
    count_radii=np.array([r2*(1+eps)**2 for eps,B in inputs]); count_bounds=np.array([B for eps,B in inputs]); lambdas=cp.Variable(len(inputs),nonneg=True) if inputs else None
    radial=np.linspace(r2,R2,nodes); scale=np.full(nodes,np.exp(-np.pi*r2/2))
    ell_scale=np.array([scale for ell in range(angular+1)])
    norms=np.exp(gammaln(np.arange(degree+1)+a)-gammaln(a)-gammaln(np.arange(degree+1)+1))
    assert seed, "Fixed radial optimization requires --seed"
    seed_data=json.load(open(seed))[0]
    coeff=cp.Constant(np.asarray(seed_data['coefficients']))
    blocks=[cp.Variable((nodes+1,nodes+1),symmetric=True)]+[cp.Variable((nodes,nodes),symmetric=True) for _ in range(angular)]
    core=[blocks[0][1:,1:]]+blocks[1:]
    vector=cp.hstack([cp.reshape(B,(nodes*nodes,),order='C') for B in core])
    cons=[B>>1e-6*np.eye(B.shape[0]) for B in blocks]
    signs=(-1.)**np.arange(degree+1)
    normalization=float(signs@coeff.value)
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
    # Radial function is frozen. Its full signs are checked separately with exact guards.
    def weights(u,ell=0):
        z=(2*np.sqrt(u)-np.sqrt(R2)-np.sqrt(r2))/(np.sqrt(R2)-np.sqrt(r2))
        return np.polynomial.chebyshev.chebvander(z,nodes-1).ravel()*scale
    def harmonic(u,v,w):
        x=(u+v-w)/(2*np.sqrt(u*v))
        q=[1.,x]
        for ell in range(1,angular):
            q.append(((2*ell+n-2)*x*q[-1]-ell*q[-2])/(ell+n-2))
        return np.asarray(q)
    def kernel_row(u,v,w):
        q=harmonic(u,v,w)
        return np.outer(q,np.outer(weights(u),weights(v)).ravel()).ravel()
    def cyclic_row(u,v,w):
        return (kernel_row(u,v,w)+kernel_row(u,w,v)+kernel_row(v,w,u))*np.exp(np.pi*(u+v+w)/3)
    def far_row(u,v,z):
        cutoff=r2 if pairwise else R2
        t=-1+z*(1+(u+v-cutoff)/(2*np.sqrt(u*v)))
        w=u+v-2*np.sqrt(u*v)*t
        return kernel_row(u,v,w)*np.exp(np.pi*(u+v)/2)
    def pair_con(u):
        w=weights(u)*np.exp(np.pi*u/2)
        diag=sum(cp.sum(cp.multiply(np.outer(weights(u,ell),weights(u,ell))*np.exp(np.pi*u),B)) for ell,B in enumerate(core))
        cross=2*np.exp(np.pi*u/2)*(w@blocks[0][1:,0])
        allowance=np.exp(np.pi*u)*((u<count_radii-1e-14).astype(float)@lambdas) if inputs else 0
        return lag(np.pi*u)@coeff+cross+diag-allowance<=-1e-6
    cons += [pair_con(u) for u in np.unique(np.r_[np.linspace(r2,R2,32*(nodes-1)+1),count_radii]) if r2<=u<=R2]
    fine=np.linspace(r2,R2,subdivisions*(nodes-1)+1)
    far_points=[(u,v,z) for i,u in enumerate(fine) for v in fine[i:] for z in np.linspace(0,1,41)]
    cyc_points=list(itertools.combinations_with_replacement(fine,3)) if not pairwise else []
    far_mat=sparse.csr_matrix(np.array([far_row(*p) for p in far_points]))
    cons += [far_mat@vector<=-1e-6]
    if cyc_points: cons += [sparse.csr_matrix(np.array([cyclic_row(*p) for p in cyc_points]))@vector<=-1e-6]
    # Polynomial Bernstein certificate on the entire antipodal far boundary.
    # Radius coordinate x=(r-rmin)/(R-rmin), where features are T_k(2x-1).
    polys=[]
    for k in range(nodes):
        polys.append(np.polynomial.Chebyshev.basis(k).convert(kind=np.polynomial.Polynomial)(np.polynomial.Polynomial([-1.,2.])))
    def antipodal_basis(lo,hi):
        out=np.zeros((nodes,nodes))
        for k,pol in enumerate(polys):
            pc=pol(np.polynomial.Polynomial([lo,hi-lo])).coef
            for i in range(nodes):
                out[i,k]=sum(pc[j]*math.comb(i,j)/math.comb(nodes-1,j) for j in range(min(i,len(pc)-1)+1))
        return out
    subdivisions_boundary=8
    bbs=[antipodal_basis(i/subdivisions_boundary,(i+1)/subdivisions_boundary)*scale for i in range(subdivisions_boundary)]
    antipodal_rows=[]
    for i,A in enumerate(bbs):
        for B in bbs[i:]:
            for arow in A:
                for brow in B:
                    base=np.outer(arow,brow).ravel()
                    antipodal_rows.append(np.outer((-1.)**np.arange(angular+1),base).ravel())
    cons += [sparse.csr_matrix(np.array(antipodal_rows))@vector <= -1e-9]
    assert n==6 and r2==1 and R2==2.25, "This local-box experiment is specific to the 6D witness"
    # Exact verifier identified a missed far lobe near equal short radii.
    # Impose a complete tensor-Bernstein sign on a surrounding valid far box.
    angular_polys=[np.polynomial.Polynomial([1.]),np.polynomial.Polynomial([-1.,16/9])]
    for ell in range(1,angular):
        angular_polys.append(((2*ell+n-2)*angular_polys[1]*angular_polys[-1]-ell*angular_polys[-2])/(ell+n-2))
    for box in [(.04,.16,.04,.16,.12,.19), (0.,.05,.82,.92,0.,.06)]:
        xl,xh,yl,yh,zl,zh=box
        ax=antipodal_basis(xl,xh)*scale;ay=antipodal_basis(yl,yh)*scale
        az=np.zeros((angular+1,angular+1))
        for ell,pol in enumerate(angular_polys):
            co=pol(np.polynomial.Polynomial([zl,zh-zl])).coef
            for i in range(angular+1):
                az[i,ell]=sum(co[j]*math.comb(i,j)/math.comb(angular,j)for j in range(min(i,len(co)-1)+1))
        local_rows=[np.outer(zrow,np.outer(xrow,yrow).ravel()).ravel()for xrow in ax for yrow in ay for zrow in az]
        cons += [sparse.csr_matrix(np.array(local_rows))@vector<=-1e-9]
    objective=cp.sum(coeff)+blocks[0][0,0]+(count_bounds@lambdas if inputs else 0)
    T=np.polynomial.Polynomial([0.,1.]); gegen=[np.polynomial.Polynomial([1.]),T]
    for ell in range(1,angular):
        gegen.append(((2*ell+n-2)*T*gegen[-1]-ell*gegen[-2])/(ell+n-2))
    start=time.time(); history=[]; saved_cuts=[];
    if seed:
        prior=json.load(open(seed))[0]
        for kind,point in prior.get('separation_cuts',[]):
            if kind=='pair':
                cons += [pair_con(float(point[0]))]
            else:
                fn=far_row if kind=='far' else cyclic_row
                cons += [fn(*point)@vector<=-2e-6]
        saved_cuts=prior.get('separation_cuts',[])
    rng=np.random.default_rng(20261111+n)
    # Closest feasible kernel to the supplied near-certificate, with a small
    # objective allowance. This avoids wandering along nearly flat SDP faces.
    proximal_terms=[]
    for ell,B in enumerate(blocks):
        sc=np.r_[1.,scale] if ell==0 else scale
        previous=np.asarray(seed_data['matrices'][ell])/np.outer(sc,sc)
        denom=max(float(np.linalg.norm(previous)),.01)
        proximal_terms.append(cp.sum_squares((B-previous)/denom))
    if inputs: proximal_terms.append(100*cp.sum_squares(lambdas-np.asarray(seed_data['count_weights'])))
    proximal=1e6*sum(proximal_terms)
    for rnd in range(rounds):
        p=cp.Problem(cp.Minimize(proximal),cons+[objective<=5.0816])
        p.solve(solver='CLARABEL',max_iter=150,tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10)
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
                pairs=[(u,v)for ii,u in enumerate(fine)for v in fine[ii:]]
                for rr in np.linspace(np.sqrt(r2),np.sqrt(R2),101):
                    uu=rr*rr;pairs.extend([(uu,uu),(r2,uu),(uu,R2)])
                for u,v in pairs:
                        outer=np.outer(weights(u),weights(v))
                        cs=[np.sum(outer*B) for ell,B in enumerate(Bvals)]
                        pol=sum((c*g for c,g in zip(cs,gegen)),np.polynomial.Polynomial([0.]))
                        cu=(u+v-(r2 if pairwise else R2))/(2*np.sqrt(u*v))
                        ts=[-1.,cu]+[float(z.real) for z in pol.deriv().roots() if abs(z.imag)<1e-7 and -1<z.real<cu]
                        extrema.extend((u,v,(t+1)/(cu+1)) for t in ts)
                starts=np.vstack([starts,np.asarray(extrema)])
                for q in extrema:
                    value=func(*q)@val
                    if value>-2e-7:
                        cons += [func(*q)@vector<=-2e-6]; saved_cuts.append([kind,list(map(float,q))]);added+=1
            values=np.asarray([func(*q)@val for q in starts])
            candidates=starts[np.argsort(values)[-32:]]
            peaks=[]
            for q in candidates:
                opt=minimize(lambda z:-1e4*(func(*z)@val),q,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-14,'gtol':1e-10,'maxiter':100})
                peaks.append((float(-opt.fun)/1e4,opt.x))
            worst[kind]=max(x[0] for x in peaks)
            for value,point in peaks:
                if value>-2e-7:
                    cons += [func(*point)@vector<=-2e-6]; saved_cuts.append([kind,list(map(float,point))]); added+=1
        def pair_value(u):
            w=weights(u)
            correction=2*w@np.asarray(blocks[0].value)[1:,0]+sum(w@np.asarray(B.value)@w for B in core)
            allowance=float((u<count_radii-1e-14).astype(float)@lambdas.value) if inputs else 0.
            return float(lag(np.pi*u)@coeff.value+np.exp(np.pi*u)*(correction-allowance))
        edges=np.unique(np.r_[np.linspace(r2,R2,64*(nodes-1)+1),count_radii])
        pair_peaks=[(pair_value(float(u)),float(u)) for u in edges]
        for lo,hi in zip(edges[:-1],edges[1:]):
            opt=minimize_scalar(lambda u:-pair_value(u),bounds=(lo,hi),method='bounded',options={'xatol':1e-14})
            pair_peaks.append((float(-opt.fun),float(opt.x)))
        worst['pair']=max(z[0] for z in pair_peaks)
        for value,u in pair_peaks:
            if value>-5e-7:
                cons += [pair_con(u)];saved_cuts.append(['pair',[u]]);added+=1
        h={'round':rnd,'objective':float(objective.value),'density_candidate':float(objective.value)/normalization*np.pi**a/math.gamma(a+1)*(r2/4)**a,'status':str(p.status),'added':added,'worst_sampled':worst,'seconds':time.time()-start}
        history.append(h); print(h,flush=True)
        result={'dimension':n,'certified':False,'warning':'Uncertified Chebyshev-radial cyclic candidate; complete exact verification is mandatory.','pairwise_control':pairwise,'radial_degree':degree,'radial_nodes':nodes,'angular_degree':angular,'exclusion_squared':str(r2),'truncation_squared':str(R2),'kernel_family':'chebyshev_radial_cyclic', 'count_squared_radii':count_radii.tolist(), 'count_bounds':count_bounds.tolist(), 'count_weights':lambdas.value.tolist() if inputs else [],'status':str(p.status),'history':history,'separation_cuts':saved_cuts}
        sc=np.r_[1,scale]
        result.update(solver_iterations=p.solver_stats.num_iters, solver_seconds=p.solver_stats.solve_time, coefficients=coeff.value.tolist(),matrices=[(blocks[0].value*np.outer(sc,sc)).tolist()]+[(B.value*np.outer(ell_scale[ell],ell_scale[ell])).tolist() for ell,B in enumerate(blocks[1:],start=1)],objective=float(objective.value),density_candidate=h['density_candidate'],max_constraint_violation=max(float(np.max(c.violation())) for c in cons if c.value() is not None))
        if output: Path(output).write_text(json.dumps([result],indent=2))
        if not added: break
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--dimension',type=int,default=17);ap.add_argument('--degree',type=int,default=31);ap.add_argument('--nodes',type=int,default=5);ap.add_argument('--angular',type=int,default=16);ap.add_argument('--ratio',type=float,default=1.2);ap.add_argument('--subdivisions',type=int,default=2);ap.add_argument('--rounds',type=int,default=5);ap.add_argument('--pairwise',action='store_true');ap.add_argument('--output',default='candidate.json')
    ap.add_argument('--exclusion',type=float,default=1.0);ap.add_argument('--counts',action='store_true');ap.add_argument('--seed');
    args=vars(ap.parse_args());args['n']=args.pop('dimension');solve(**args)
