"""Exact rational construction; outward-rounded Bernstein exhaustive covering.
No sampled point is used for proof acceptance. Frozen inputs are not modified.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,time,math,sys
from fractions import Fraction
import numpy as np
from flint import fmpq,fmpq_poly,fmpq_mat,fmpq_mpoly_ctx,arb,ctx
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'density_kissing16_reweight_frozen_20261111';ctx.prec=192
SHAS={'candidate.json':'b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573','factors.json':'6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954'}
for name,sha in SHAS.items():assert hashlib.sha256((SRC/name).read_bytes()).hexdigest()==sha
q=lambda x:fmpq(str(Fraction(str(x))))
DATA=json.loads((SRC/'factors.json').read_text(),parse_float=str);W=json.loads((SRC/'candidate.json').read_text(),parse_float=str)['record'];TH=list(map(q,W['theta']));RR=list(map(q,W['local_angular_coefficients']));CAP=fmpq(5003,10000)
MC=fmpq_mpoly_ctx.get(['u','v','t']);u,v,t=MC.gens();start=time.monotonic()
def gegen(x,n,d):
    p=[x*0+1]
    if d:p.append(x)
    for k in range(1,d):p.append(((2*k+n-2)*x*p[-1]-k*p[-2])/fmpq(k+n-2))
    return p

def init():
    FF=[]
    for block in DATA['blocks']:
        n=block['size'];L=fmpq_mat(n,n)
        for i,j,s in block['L_entries']:L[i,j]=q(s)
        F=L.transpose()*L
        for i in range(n):F[i,i]+=q(block['lambda_decimal'])
        FF.append(F)
    BB=FF[16];FF=FF[17:]
    g=t-u*v;z=(1-u*u)*(1-v*v);Q=[MC.constant(1),g]
    for k in range(1,16):Q.append(((2*k+13)*g*Q[-1]-k*z*Q[-2])/fmpq(k+13))
    X=fmpq_poly([0,1]);ker=MC.constant(BB[1,1])
    for k,F in enumerate(FF):
        if not TH[k]:continue
        n=17-k;P=gegen(X,16+2*k,16-k);M=fmpq_mat([[P[i][j] for i in range(n)] for j in range(n)]);Z=M*F*M.transpose()
        bp=MC.from_dict({(i,j,0):Z[i,j] for i in range(n) for j in range(n) if Z[i,j]})
        ker+=TH[k]*bp*Q[k]
    poly=ker+ker.compose(u,t,v)+ker.compose(v,t,u)
    Pu=gegen(u,16,24)
    rpoly=sum((RR[k-1]*Pu[k] for k in range(1,25)),MC.constant(0))
    poly+=rpoly+rpoly.compose(v,u,t)+rpoly.compose(t,v,u)
    # Exact rational verification of permutation symmetry and direct recurrence at distinct points.
    assert poly==poly.compose(v,u,t)==poly.compose(u,t,v)
    for point in [(fmpq(-3,5),fmpq(1,7),fmpq(1,3)),(fmpq(-1),-CAP,CAP),(CAP,CAP,CAP)]:
        def kernel(a,b,c):
            gg=c-a*b;zz=(1-a*a)*(1-b*b);qq=[fmpq(1),gg]
            for k in range(1,16):qq.append(((2*k+13)*gg*qq[-1]-k*zz*qq[-2])/fmpq(k+13))
            val=BB[1,1]
            for k,F in enumerate(FF):
                if not TH[k]:continue
                pa=gegen(a,16+2*k,16-k);pb=gegen(b,16+2*k,16-k)
                val+=TH[k]*qq[k]*sum((pa[i]*F[i,j]*pb[j] for i in range(len(pa)) for j in range(len(pb))),fmpq(0))
            return val
        a,b,c=point;value=kernel(a,b,c)+kernel(a,c,b)+kernel(b,c,a)
        for w in point:
            ps=gegen(w,16,24);value+=sum((RR[k-1]*ps[k] for k in range(1,25)),fmpq(0))
        assert poly(*point)==value
    D=poly.to_dict();d=int(max(max(e) for e in D));n=d+1;print('EXACT_POLYNOMIAL',len(D),'axisdegree',d,'totaldegree',poly.total_degree(),'seconds',time.monotonic()-start,flush=True)
    (ROOT/'polynomial.json').write_text(json.dumps({'source_sha256':SHAS,'terms':[[list(map(int,e)),str(z)]for e,z in D.items()]},separators=(',',':')))
    # Substitute x=-1+(1+cap)s, then power-to-degree-d Bernstein.
    a=fmpq(-1);delta=CAP+1
    BT=fmpq_mat([[sum((fmpq(math.comb(e,j)*math.comb(i,j),math.comb(d,j))*a**(e-j)*delta**j for j in range(min(i,e)+1)),fmpq(0)) for e in range(n)] for i in range(n)])
    A=np.empty((n,n,n),dtype=object);A.fill(fmpq(0))
    for e,z in D.items():A[e]=z
    for axis in range(3):
        V=np.moveaxis(A,axis,0).reshape(n,n*n);M=BT*fmpq_mat(V.tolist());V=np.array([[M[i,j] for j in range(n*n)]for i in range(n)],dtype=object);A=np.moveaxis(V.reshape(n,n,n),0,axis)
        print('TRANSFORM',axis,'seconds',time.monotonic()-start,flush=True)
    lo=np.empty(A.shape);hi=np.empty(A.shape)
    for e in np.ndindex(A.shape):
        av=arb(A[e]);lo[e]=np.nextafter(float(av.lower()),-np.inf);hi[e]=np.nextafter(float(av.upper()),np.inf)
    np.savez_compressed(ROOT/'bernstein_initial.npz',lo=lo,hi=hi)
    print('INITIAL_BOUNDS',float(np.min(lo)),float(np.max(hi)),'maxwidth',float(np.max(hi-lo)),flush=True)
    return lo,hi

def split(l,h,axis):
    a=np.moveaxis(l,axis,0).copy();b=np.moveaxis(h,axis,0).copy();d=a.shape[0]-1
    ll=np.empty_like(a);lh=np.empty_like(b);rl=np.empty_like(a);rh=np.empty_like(b)
    ll[0]=a[0];lh[0]=b[0];rl[d]=a[-1];rh[d]=b[-1]
    for k in range(1,d+1):
        a=np.nextafter(np.nextafter(a[:-1]+a[1:],-np.inf)*.5,-np.inf);b=np.nextafter(np.nextafter(b[:-1]+b[1:],np.inf)*.5,np.inf)
        ll[k]=a[0];lh[k]=b[0];rl[d-k]=a[-1];rh[d-k]=b[-1]
    return tuple(np.moveaxis(z,0,axis)for z in (ll,lh,rl,rh))

def gram_exclude(box):
    # Exact rational range enclosure for 1+2uvt-u²-v²-t².
    # Only invoked for boxes whose cyclic upper bound is nonnegative.
    ints=[(Fraction(-1)+Fraction(15003,10000)*a,Fraction(-1)+Fraction(15003,10000)*b) for a,b in box]
    product=max(a*b*c for a in ints[0] for b in ints[1] for c in ints[2])
    min_squares=sum(0 if a<=0<=b else min(a*a,b*b) for a,b in ints)
    return 1+2*product-min_squares<0

def run(limit=50000,seconds=600):
    if '--reuse' in sys.argv:
        a=np.load(ROOT/'bernstein_initial.npz');lo,hi=a['lo'],a['hi']
    else:lo,hi=init()
    stack=[(lo,hi,((Fraction(0),Fraction(1)),)*3,0)];checked=passed=excluded_sym=excluded_gram=0;maxdepth=0;fail=None;peakstack=1
    while stack:
        l,h,box,depth=stack.pop();checked+=1;maxdepth=max(maxdepth,depth)
        if np.max(h)<0:passed+=1;continue
        if box[0][0]>box[1][1] or box[1][0]>box[2][1]:excluded_sym+=1;continue
        if gram_exclude(box):excluded_gram+=1;continue
        if checked>=limit or depth>=100 or time.monotonic()-start>seconds:
            fail={'box':[[str(a),str(b)]for a,b in box],'depth':depth,'upper':float(np.max(h)),'remaining':len(stack),'reason':'bounded work limit'};break
        mid=(l+h)*.5;scores=[float(np.max(np.abs(np.diff(mid,axis=i)))) for i in range(3)]
        # Balance physical widths to prevent pathological near-flat splits.
        axis=int(np.argmax(scores));ll,lh,rl,rh=split(l,h,axis);a,b=box[axis];m=(a+b)/2;bx1=list(box);bx2=list(box);bx1[axis]=(a,m);bx2[axis]=(m,b)
        stack.extend([(rl,rh,tuple(bx2),depth+1),(ll,lh,tuple(bx1),depth+1)]);peakstack=max(peakstack,len(stack))
        if checked%1000==0:print('PROGRESS',checked,'passed',passed,'gram',excluded_gram,'sym',excluded_sym,'stack',len(stack),'depth',maxdepth,'seconds',time.monotonic()-start,flush=True)
    out={'verified':fail is None and not stack,'source_sha256':SHAS,'checked':checked,'passed':passed,'excluded_symmetry':excluded_sym,'excluded_gram':excluded_gram,'max_depth':maxdepth,'peak_stack':peakstack,'seconds':time.monotonic()-start,'failure':fail,'method':'Exact rational cyclic polynomial and exact rational tensor-Bernstein initialization on [-1,5003/10000]^3. Outward double dyadic de Casteljau subdivision. Exact rational Gram determinant interval exclusion. Exact polynomial permutation symmetry; restrict to u<=v<=t. No samples used for acceptance.'}
    (ROOT/'bernstein_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 and sys.argv[1].isdigit() else 50000,600)
