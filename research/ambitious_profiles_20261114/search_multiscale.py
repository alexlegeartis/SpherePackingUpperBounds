"""Larger matrix-valued cyclic SDP with three or more radial Gaussian features.

The long-edge cap is controlled by exact Bernstein coefficients in the
common fast/slow ratio. This is a sampled discovery problem, never a proof.
"""
import argparse,json,math,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np,cvxpy as cp
from scipy.special import eval_genlaguerre,gammaln
from flint import arb,arb_poly,ctx
from common import HERE,source,base_angular

ap=argparse.ArgumentParser();ap.add_argument('dimension',type=int)
ap.add_argument('--rates',default='40,50,60');ap.add_argument('--degree',type=int,default=63)
ap.add_argument('--s',default='6');ap.add_argument('--cut',type=float,default=40);ap.add_argument('--round',type=int,default=0)
ap.add_argument('--fixed-radial')
args=ap.parse_args();n,d=args.dimension,args.degree;s=F(args.s)
rates=list(map(F,args.rates.split(',')));m=len(rates)
assert m>=2 and rates[0]>=40 and all(x<y for x,y in zip(rates,rates[1:]))
gamma=rates[1]-rates[0];powers=[int((b-rates[0])/gamma) for b in rates]
assert all(F(k)==(b-rates[0])/gamma for k,b in zip(powers,rates))
K=max(powers);D=HERE/f'multiscale_n{n}_m{m}_d{d}_s{str(s).replace("/","_")}_b{str(args.rates).replace(",","_")}'
if args.fixed_radial:D=Path(str(D)+'_fr')
D.mkdir(exist_ok=True);H=D/f'round{args.round:02}';H.mkdir(exist_ok=True)
(H/'search_source.py').write_text(Path(__file__).read_text())
ctx.prec=512;start=time.time()
refpath=base_angular(n);ref=json.loads(refpath.read_text())
bell=np.array([float(F(x)) for x in ref['gegenbauer_coefficients']]);B=float(F(ref['B_exact']))
weights=bell/B;ellmax=len(bell)-1;active=[ell for ell,z in enumerate(bell) if z]
p,h,ls=source(n,d);AP=arb_poly([arb(z) for z in p.coeffs()]);AH=arb_poly([arb(z) for z in h.coeffs()])
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
indices=[(i,j) for i in range(m) for j in range(i,m)];pairs=len(indices)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*np.asarray(x))/norm[k] for k in range(d+1)]).T
def pol(poly,x):return np.array([float(poly(arb(float(z))).mid()) for z in x])
def gegen(t):
    t=np.asarray(t);g=[np.ones_like(t),t]
    for ell in range(1,ellmax):g.append(((2*ell+n-2)*t*g[-1]-ell*g[-2])/(ell+n-2))
    return np.stack(g,axis=-1)
def tri_rows(t,u,T):
    t,u,T=np.broadcast_arrays(np.asarray(t),np.asarray(u),np.asarray(T))
    w2=1+u*u-2*t*u;w=np.sqrt(w2)
    g=[gegen(t),gegen((1-u*t)/w),gegen((u-t)/w)]
    radius=[T,T*u*u,T*w2]
    rat=[np.stack([np.exp(-float(b-rates[0])*(v-1)) for b in rates],axis=-1) for v in radius]
    gaussian=[np.ones_like(t),np.exp(-float(rates[0])*T*(w2-u*u)),np.exp(-float(rates[0])*T*(w2-1))]
    rows=np.zeros((len(t),len(active),pairs))
    for vertex,(left,right) in enumerate([(0,1),(0,2),(1,2)]):
        for ij,(i,j) in enumerate(indices):
            factor=rat[left][:,i]*rat[right][:,j]
            if i!=j:factor+=rat[left][:,j]*rat[right][:,i]
            rows[:,:,ij]+=(gaussian[vertex]*factor)[:,None]*g[vertex][:,active]
    return rows.reshape(len(t),-1)*np.repeat(weights[active],pairs)
def cap_rows(t):
    gv=gegen(t)[:,active]*weights[active]
    rows=[]
    for j in range(K+1):
        coeff=np.zeros(pairs);coeff[indices.index((0,0))]=1
        for k in range(1,m):
            if j>=powers[k]:coeff[indices.index((0,k))]=math.comb(j,powers[k])/math.comb(K,powers[k])
        rows.append((gv[:,:,None]*coeff).reshape(len(t),-1))
    return np.vstack(rows)
tt,uu=np.meshgrid(np.r_[np.linspace(-1,0,31),np.linspace(0,.5,51)],np.r_[np.linspace(1,1.3,21),np.linspace(1.3,2,13),np.linspace(2,5,9)])
ts,us=tt.ravel(),uu.ravel();mask=2*ts*us<=1;ts,us=ts[mask],us[mask]
tri=np.vstack([tri_rows(ts,us,np.full_like(ts,T)) for T in [1,1.025,1.075,1.15,1.3,1.6,2,3]])
cap=cap_rows(np.linspace(-1,.5001,2001))
cutsfile=D/'cuts.json';cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {}
if cuts.get('triangles'):
    z=np.array(cuts['triangles']);tri=np.vstack([tri,tri_rows(z[:,0],z[:,1],z[:,2])])
if cuts.get('cap'):cap=np.vstack([cap,cap_rows(np.array(cuts['cap']))])
for rows in [tri,cap]:rows/=np.maximum(abs(rows).max(axis=1,keepdims=True),1e-200)
X=cp.Variable((m+1,m+1),symmetric=True);A={0:X[1:,1:]}
for ell in active[1:]:A[ell]=cp.Variable((m,m),symmetric=True)
kernel=cp.hstack([cp.hstack([A[ell][i,j] for i,j in indices]) for ell in active])
if args.fixed_radial:
    fixed=json.loads(Path(args.fixed_radial).read_text())
    fixed_polynomial='physical_polynomial' in fixed
    if fixed_polynomial:
        fixed_H0=F(fixed['fourier_polynomial'][0]);assert fixed_H0>0
        fixed_coeffs=[F(z)/fixed_H0 for z in fixed['physical_polynomial']]
        fixed_AP=arb_poly([arb(z.numerator)/z.denominator for z in fixed_coeffs])
        fixed_f0=float(fixed_AP(arb(0)).mid());fixed_values=[0.]*(d+2)
    else:fixed_values=fixed.get('radial',fixed.get('x',[])[:-1])
    assert fixed['dimension']==n and len(fixed_values)<=d+2
    fixed_values=list(fixed_values)+[0.]*(d+2-len(fixed_values))
    radial=cp.Constant(np.array(fixed_values))
else:radial=cp.Variable(d+2)
cons=[X>>0,X[0,1]<=-1e-5,A[0][0,0]>=1e-5,tri@kernel<=-1e-7,cap@kernel<=-1e-7]
if not args.fixed_radial:cons.append(radial[0]>=0)
for ell in active[1:]:cons.extend([A[ell]>>0,A[ell][0,0]==A[0][0,0]])
limit=max(300.,4*d+n+20.)
fouriergrid=np.r_[np.linspace(0,limit,4001),np.geomspace(limit,2500,301),cuts.get('fourier',[])]
physicalgrid=np.r_[np.linspace(args.cut,limit,3001),np.geomspace(limit,2500,301),cuts.get('physical',[])]
pairgrid=np.unique(np.r_[np.linspace(float(s)*np.pi,args.cut,1201),float(s)*np.pi+np.geomspace(1e-7,2.,151),cuts.get('pair',[])])
fr=np.c_[-pol(AH,fouriergrid),-lag(fouriergrid)*sign]
pr=np.c_[pol(AP,physicalgrid),lag(physicalgrid)]
fr/=np.maximum(abs(fr).max(axis=1,keepdims=True),1e-200);pr/=np.maximum(abs(pr).max(axis=1,keepdims=True),1e-200)
if not args.fixed_radial:
    cons.extend([fr@radial<=-1e-8,pr@radial<=-1e-8,np.r_[float(h[0]),sign]@radial==1])
v=pairgrid/(float(s)*np.pi);hs=np.stack([np.exp(-float(b)*(v-1)) for b in rates],axis=1)
phys=np.c_[pol(AP,pairgrid),lag(pairgrid)]*(B*np.exp(-pairgrid))[:,None]
hd=np.stack([(1 if i==j else 2)*hs[:,i]*hs[:,j] for i,j in indices],axis=1)
diag=np.tile(hd,(1,len(active)))*np.repeat(weights[active],pairs)
scale=np.maximum(np.max(abs(phys),axis=1),np.maximum(np.max(abs(diag),axis=1),2*np.max(hs,axis=1)))
if args.fixed_radial and fixed_polynomial:
    fixed_pair=pol(fixed_AP,pairgrid)*B*np.exp(-pairgrid)
    scale=np.maximum(abs(fixed_pair),np.maximum(np.max(abs(diag),axis=1),2*np.max(hs,axis=1)))
    physical_expression=fixed_pair/scale
else:physical_expression=(phys/scale[:,None])@radial
cons.append(physical_expression+(diag/scale[:,None])@kernel+2*(hs/scale[:,None])@X[0,1:]<=-1e-9)
objective=(fixed_f0 if args.fixed_radial and fixed_polynomial else np.r_[float(p[0]),np.ones(d+1)]@radial)+X[0,0]
problem=cp.Problem(cp.Minimize(objective),cons);attempts=[]
for tolerance in [1e-10,1e-8,1e-7]:
    try:
        value=problem.solve(solver='CLARABEL',tol_feas=tolerance,tol_gap_abs=tolerance,tol_gap_rel=tolerance,
                            max_iter=200,max_threads=2,verbose=False)
        attempts.append({'tolerance':tolerance,'status':problem.status})
        if X.value is not None:break
    except cp.error.SolverError as exc:attempts.append({'tolerance':tolerance,'error':str(exc)})
(H/'solver_attempts.json').write_text(json.dumps(attempts,indent=2)+'\n')
assert X.value is not None and radial.value is not None,attempts
blocks=[np.zeros((m,m)).tolist() for _ in range(ellmax+1)]
for ell in active:blocks[ell]=A[ell].value.tolist()
out={'dimension':n,'features':m,'degree':d,'s_exact':str(s),'rates_exact':list(map(str,rates)),
     'powers':powers,'fixed_slow_angular':True,'cut_X_exact':str(F(str(args.cut))),
     'reference_angular':str(refpath.relative_to(HERE.parents[1])),
     'fixed_radial_input':args.fixed_radial,
     'radial':radial.value.tolist(),'origin_block':X.value.tolist(),'blocks':blocks,
     'objective':float(value),'density_candidate':float(value)*np.pi**(n/2)*(float(s)/4)**(n/2)/math.gamma(n/2+1),
     'certified':False,'solver_status':problem.status,'solver_attempts':attempts,'seconds':time.time()-start,
     'minimum_eigenvalue':float(min(np.linalg.eigvalsh(X.value).min(),*(np.linalg.eigvalsh(A[ell].value).min() for ell in active[1:]))),
     'max_triangle_residual':float(max(tri@kernel.value)),'max_cap_residual':float(max(cap@kernel.value)),
     'origin_cost':float(X.value[0,0]),'cross_vector':X.value[0,1:].tolist(),
     'triangle_rows':len(tri),'cap_rows':len(cap),'pair_rows':len(pairgrid)}
(H/'candidate.json').write_text(json.dumps(out,indent=2)+'\n');(D/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ['radial','origin_block','blocks']},flush=True)
