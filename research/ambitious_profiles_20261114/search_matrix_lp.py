"""PSD-kernel discovery by nonnegative Gram rays and LP column generation.

Every cone column is a rank-one PSD matrix. Negative eigenvectors of sampled
dual slack matrices add new columns. All outputs remain numerical research.
"""
import argparse,json,math,time,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
from scipy.sparse import csr_matrix,hstack,vstack,block_diag
from flint import arb,arb_poly,ctx
from common import HERE,source,base_angular

ap=argparse.ArgumentParser();ap.add_argument('dimension',type=int)
ap.add_argument('--rates',default='40,60');ap.add_argument('--degree',type=int,default=63)
ap.add_argument('--s',default='6');ap.add_argument('--cut',type=float,default=40)
ap.add_argument('--iterations',type=int,default=15);ap.add_argument('--round',type=int,default=0)
ap.add_argument('--cuts-from');ap.add_argument('--fixed-polynomial')
args=ap.parse_args();n,d=args.dimension,args.degree;s=F(args.s)
rates=list(map(F,args.rates.split(',')));m=len(rates);assert m>=2 and rates[0]>=40
assert all(x<y for x,y in zip(rates,rates[1:]))
gamma=rates[1]-rates[0];powers=[int((b-rates[0])/gamma) for b in rates]
assert all(F(k)==(b-rates[0])/gamma for k,b in zip(powers,rates))
K=max(powers)
D=HERE/f'rays_n{n}_m{m}_d{d}_s{str(s).replace("/","_")}_b{args.rates.replace(",","_")}'
if args.fixed_polynomial:D=Path(str(D)+'_fp')
D.mkdir(exist_ok=True);H=D/f'round{args.round:02}';H.mkdir(exist_ok=True)
(H/'search_source.py').write_text(Path(__file__).read_text());ctx.prec=512;start=time.time()
refpath=base_angular(n);ref=json.loads(refpath.read_text())
bell=np.array([float(F(x)) for x in ref['gegenbauer_coefficients']]);B=float(F(ref['B_exact']))
weights=bell/B;ellmax=len(bell)-1;active=[ell for ell,z in enumerate(bell) if z]
indices=[(i,j) for i in range(m) for j in range(i,m)];q=len(indices)
origin_indices=[(i,j) for i in range(m+1) for j in range(i,m+1)];q0=len(origin_indices)
dimensions=[m+1]+[m]*(len(active)-1)
slices=[slice(0,q0)]+[slice(q0+(k-1)*q,q0+k*q) for k in range(1,len(active))]
N=q0+(len(active)-1)*q
mapping=np.zeros((len(active)*q,N))
for k,ell in enumerate(active):
    for ij,(i,j) in enumerate(indices):
        col=origin_indices.index((i+1,j+1)) if k==0 else q0+(k-1)*q+ij
        mapping[k*q+ij,col]=1
mapping=csr_matrix(mapping)
def into_full(x):return np.asarray(x@mapping)
def normalize(x):return x/np.maximum(abs(x).max(axis=1,keepdims=True),1e-200)
def gegen(t):
    t=np.asarray(t);g=[np.ones_like(t),t]
    for ell in range(1,ellmax):g.append(((2*ell+n-2)*t*g[-1]-ell*g[-2])/(ell+n-2))
    return np.stack(g,axis=-1)
def tri_rows(t,u,T):
    t,u,T=np.broadcast_arrays(np.asarray(t),np.asarray(u),np.asarray(T))
    w2=1+u*u-2*t*u;w=np.sqrt(w2)
    gs=[gegen(t),gegen((1-u*t)/w),gegen((u-t)/w)]
    radii=[T,T*u*u,T*w2]
    rr=[np.stack([np.exp(-float(b-rates[0])*(v-1)) for b in rates],axis=-1) for v in radii]
    ww=[np.ones_like(t),np.exp(-float(rates[0])*T*(w2-u*u)),np.exp(-float(rates[0])*T*(w2-1))]
    rows=np.zeros((len(t),len(active),q))
    for vertex,(left,right) in enumerate([(0,1),(0,2),(1,2)]):
        for ij,(i,j) in enumerate(indices):
            factor=rr[left][:,i]*rr[right][:,j]
            if i!=j:factor+=rr[left][:,j]*rr[right][:,i]
            rows[:,:,ij]+=(ww[vertex]*factor)[:,None]*gs[vertex][:,active]
    return into_full(rows.reshape(len(t),-1)*np.repeat(weights[active],q))
def cap_rows(t):
    gv=gegen(t)[:,active]*weights[active];rows=[]
    for j in range(K+1):
        cf=np.zeros(q);cf[indices.index((0,0))]=1
        for i in range(1,m):
            if j>=powers[i]:cf[indices.index((0,i))]=math.comb(j,powers[i])/math.comb(K,powers[i])
        rows.append(into_full((gv[:,:,None]*cf).reshape(len(t),-1)))
    return np.vstack(rows)
cuts={}
for cutfile in ([Path(args.cuts_from)] if args.cuts_from else [])+[D/'cuts.json']:
    if cutfile.exists():
        for key,vals in json.loads(cutfile.read_text()).items():cuts.setdefault(key,[]).extend(vals)
tt,uu=np.meshgrid(np.r_[np.linspace(-1,0,13),np.linspace(0,.5,21)],
                 np.array([1.,1.025,1.05,1.1,1.2,1.4,1.8,2.5,3.5,5.]))
ts,us=tt.ravel(),uu.ravel();mask=2*ts*us<=1;ts,us=ts[mask],us[mask]
tri=np.vstack([tri_rows(ts,us,np.full_like(ts,T)) for T in [1,1.05,1.15,1.4,2,3]])
if cuts.get('triangles'):
    zz=np.array(cuts['triangles']);tri=np.vstack([tri,tri_rows(zz[:,0],zz[:,1],zz[:,2])])
cap=cap_rows(np.r_[np.linspace(-1,.5001,301),cuts.get('cap',[])])
tri,cap=normalize(tri),normalize(cap)
E=np.zeros((len(active)-1,N))
for k in range(1,len(active)):
    E[k-1,slices[k].start+indices.index((0,0))]=1
    E[k-1,origin_indices.index((1,1))]=-1
direct=np.zeros((2,N));direct[0,origin_indices.index((0,1))]=1;direct[1,origin_indices.index((1,1))]=-1
p,h,ls=source(n,d);AP=arb_poly([arb(z) for z in p.coeffs()]);AH=arb_poly([arb(z) for z in h.coeffs()])
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*np.asarray(x))/norm[k] for k in range(d+1)]).T
def pol(poly,x):return np.array([float(poly(arb(float(z))).mid()) for z in x])
pairgrid=np.unique(np.r_[np.linspace(float(s)*np.pi,args.cut,801),float(s)*np.pi+np.geomspace(1e-8,2,101),cuts.get('pair',[])])
v=pairgrid/(float(s)*np.pi);hs=np.stack([np.exp(-float(b)*(v-1)) for b in rates],axis=1)
hd=np.stack([(1 if i==j else 2)*hs[:,i]*hs[:,j] for i,j in indices],axis=1)
pairK=into_full(np.tile(hd,(1,len(active)))*np.repeat(weights[active],q))
for i in range(m):pairK[:,origin_indices.index((0,i+1))]=2*hs[:,i]
if args.fixed_polynomial:
    fixed=json.loads(Path(args.fixed_polynomial).read_text());H0=F(fixed['fourier_polynomial'][0]);assert H0>0
    cs=[F(x)/H0 for x in fixed['physical_polynomial']]
    fp=arb_poly([arb(x.numerator)/x.denominator for x in cs])
    fixedpair=pol(fp,pairgrid)*B*np.exp(-pairgrid);fixedorigin=float(cs[0]);Pdim=0
    scale=np.maximum(abs(fixedpair),np.max(abs(pairK),axis=1));pairK/=scale[:,None]
    pair_rhs=-fixedpair/scale-1e-9
    Aphys=np.zeros((len(tri)+len(cap)+2+len(pairK),0));eqphys=np.zeros((len(E),0))
    KernelI=np.vstack([tri,cap,direct,pairK]);rhs=np.r_[np.full(len(tri)+len(cap),-1e-7),[-1e-5,-1e-5],pair_rhs]
    EqK=E;eqrhs=np.zeros(len(E));cphys=np.zeros(0)
else:
    Pdim=d+2;fixedorigin=0.
    limit=max(300.,4*d+n+20.)
    yg=np.r_[np.linspace(0,limit,2001),np.geomspace(limit,1500,101),cuts.get('fourier',[])]
    xg=np.r_[np.linspace(args.cut,limit,1501),np.geomspace(limit,1500,101),cuts.get('physical',[])]
    frows=normalize(np.c_[-pol(AH,yg),-lag(yg)*sign]);prows=normalize(np.c_[pol(AP,xg),lag(xg)])
    pairP=np.c_[pol(AP,pairgrid),lag(pairgrid)]*(B*np.exp(-pairgrid))[:,None]
    scale=np.maximum(np.max(abs(pairP),axis=1),np.max(abs(pairK),axis=1));pairP/=scale[:,None];pairK/=scale[:,None]
    KernelI=np.vstack([tri,cap,direct,pairK,np.zeros((len(frows)+len(prows),N))])
    Aphys=np.vstack([np.zeros((len(tri)+len(cap)+2,Pdim)),pairP,frows,prows])
    rhs=np.r_[np.full(len(tri)+len(cap),-1e-7),[-1e-5,-1e-5],np.full(len(pairK),-1e-9),np.full(len(frows)+len(prows),-1e-8)]
    EqK=np.vstack([E,np.zeros((1,N))]);eqphys=np.vstack([np.zeros((len(E),Pdim)),np.r_[float(h[0]),sign]])
    eqrhs=np.r_[np.zeros(len(E)),1.];cphys=np.r_[float(p[0]),np.ones(d+1)]
cK=np.zeros(N);cK[origin_indices.index((0,0))]=1
KernelI=csr_matrix(KernelI);Aphys=csr_matrix(Aphys);EqK=csr_matrix(EqK);eqphys=csr_matrix(eqphys)
directions=[]
for block,size in enumerate(dimensions):
    vv=[]
    for i in range(size):
        z=np.zeros(size);z[i]=1;vv.append(z)
    for i in range(size):
        for j in range(i+1,size):
            for sg in [-1,1]:
                z=np.zeros(size);z[i]=1;z[j]=sg;vv.append(z)
    if block==0:
        for k in range(1,m):
            for t in [F(1,2),F(1),F(2)]:
                z=np.zeros(size);z[0]=1;z[1]=-1-float(t);z[k+1]=float(t);vv.append(z)
    directions.append(vv)
saved=D/'directions.json'
if saved.exists():directions=[[np.array([float(F(z)) for z in ray]) for ray in block] for block in json.loads(saved.read_text())]
def ray_map():
    blocks=[]
    for k,vv in enumerate(directions):
        ij=origin_indices if k==0 else indices
        blocks.append(np.array([[v[i]*v[j] for v in vv] for i,j in ij]))
    return block_diag(blocks,format='csr')
history=[]
for iteration in range(args.iterations):
    rays=ray_map();Ar=KernelI@rays;Er=EqK@rays
    A=hstack([Aphys,Ar],format='csr');Eq=hstack([eqphys,Er],format='csr')
    c=np.r_[cphys,np.asarray(cK@rays).ravel()]
    bounds=([(0,None)]+[(None,None)]*(Pdim-1) if Pdim else [])+[(0,None)]*rays.shape[1]
    for method in ['highs','highs-ipm']:
        result=linprog(c,A_ub=A,b_ub=rhs,A_eq=Eq,b_eq=eqrhs,bounds=bounds,method=method,
                      options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'time_limit':90})
        if result.success:break
    if not result.success:raise RuntimeError(result.message)
    flat=np.asarray(rays@result.x[Pdim:]).ravel()
    gradient=cK-KernelI.T@result.ineqlin.marginals-EqK.T@result.eqlin.marginals
    additions=0;minimum=0.
    for k,size in enumerate(dimensions):
        ij=origin_indices if k==0 else indices
        rm=np.zeros((size,size))
        for z,(i,j) in zip(gradient[slices[k]],ij):rm[i,j]=rm[j,i]=z if i==j else z/2
        eig,vec=np.linalg.eigh(rm);minimum=min(minimum,float(eig[0]))
        if eig[0]<-1e-7:
            v=vec[:,0];v=v if v[np.argmax(abs(v))]>0 else -v
            v=np.array([float(F(float(z)).limit_denominator(10**7)) for z in v])
            if not any(np.linalg.norm(v-old)<1e-7 or np.linalg.norm(v+old)<1e-7 for old in directions[k]):
                directions[k].append(v);additions+=1
    objective=float(result.fun+fixedorigin)
    row={'iteration':iteration,'objective':objective,'density_candidate':objective*np.pi**(n/2)*(float(s)/4)**(n/2)/math.gamma(n/2+1),
         'rays':rays.shape[1],'added_rays':additions,'minimum_sampled_dual_eigenvalue':minimum,'seconds':time.time()-start}
    history.append(row);print(row,flush=True)
    X=np.zeros((m+1,m+1))
    for z,(i,j) in zip(flat[slices[0]],origin_indices):X[i,j]=X[j,i]=z
    matrices=[np.zeros((m,m)).tolist() for _ in range(ellmax+1)];matrices[0]=X[1:,1:].tolist()
    for k,ell in enumerate(active[1:],start=1):
        mat=np.zeros((m,m))
        for z,(i,j) in zip(flat[slices[k]],indices):mat[i,j]=mat[j,i]=z
        matrices[ell]=mat.tolist()
    out={'dimension':n,'features':m,'degree':d,'s_exact':str(s),'rates_exact':list(map(str,rates)),
         'powers':powers,'fixed_slow_angular':True,'cut_X_exact':str(F(str(args.cut))),
         'reference_angular':str(refpath.relative_to(HERE.parents[1])),
         'fixed_radial_input':args.fixed_polynomial,'radial':result.x[:Pdim].tolist() if Pdim else [],
         'origin_block':X.tolist(),'blocks':matrices,'objective':objective,'density_candidate':row['density_candidate'],
         'certified':False,'solver_status':result.message,'solver_family':'PSD Gram-ray column generation',
         'origin_cost':float(X[0,0]),'seconds':time.time()-start,'iteration':iteration,
         'rays':rays.shape[1],'max_normalized_violation':float(max(A@result.x-rhs))}
    (H/f'iteration{iteration:02}.json').write_text(json.dumps(out,indent=2)+'\n')
    (D/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
    (H/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    # Save the expanded dictionary for a later solve; checkpoints above retain the actual matrices.
    saved.write_text(json.dumps([[[str(F(float(z))) for z in ray] for ray in block] for block in directions],indent=2)+'\n')
    if additions==0:break
