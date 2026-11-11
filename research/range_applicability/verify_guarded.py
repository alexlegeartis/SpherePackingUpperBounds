"""Fail-closed exact/interval verifier for a solid-harmonic cyclic witness.

Uses exact rational PSD, Fourier identities and Sturm tests, exact rational
Bernstein initialization, outward-rounded IEEE binary64 de Casteljau,
and Arb interval bounds.  It does not use numerical optimization.
"""
import argparse, itertools, json, math, time
from pathlib import Path
from reciprocal_guards import verify_guards
from fractions import Fraction
import numpy as np
from flint import fmpq as Q, fmpq_poly as Poly, fmpq_mat as Mat, arb, ctx
assert __debug__, 'Assertions must remain enabled.'
ctx.prec=256


def add(A,B,scale=Q(1)):
    C=A.copy()
    for e,v in B.items():
        C[e]=C.get(e,Q(0))+scale*v
        if not C[e]: del C[e]
    return C


def shift(A,e,c=Q(1)):
    return {tuple(i+j for i,j in zip(k,e)):c*v for k,v in A.items() if c*v}


def mul(A,B):
    C={}
    for e,v in B.items(): C=add(C,shift(A,e,v))
    return C


def permute(A,p): return {tuple(e[j] for j in p):v for e,v in A.items()}


def ev(p,x):
    y=arb(0)
    for c in reversed(p.coeffs()): y=y*x+arb(c)
    return y


def mag(p,x,h):
    co=[arb(0)]*max(1,p.degree()+1)
    for c in reversed(p.coeffs()):
        for j in range(p.degree(),0,-1): co[j]=co[j]*x+co[j-1]
        co[0]=co[0]*x+arb(c)
    z=arb(0)
    for c in reversed(co): z=z*h+abs(c)
    return z


def split(l,h,axis):
    a=np.moveaxis(l,axis,0).copy();b=np.moveaxis(h,axis,0).copy();d=a.shape[0]-1
    ll=np.empty_like(a);lh=np.empty_like(b);rl=np.empty_like(a);rh=np.empty_like(b)
    ll[0]=a[0];lh[0]=b[0];rl[d]=a[-1];rh[d]=b[-1]
    for k in range(1,d+1):
        a=np.nextafter(np.nextafter(a[:-1]+a[1:],-np.inf)*.5,-np.inf)
        b=np.nextafter(np.nextafter(b[:-1]+b[1:],np.inf)*.5,np.inf)
        ll[k]=a[0];lh[k]=b[0];rl[d-k]=a[-1];rh[d-k]=b[-1]
    return tuple(np.moveaxis(z,0,axis) for z in [ll,lh,rl,rh])


def bernstein(D,box):
    degrees=[max(e[k] for e in D) for k in range(3)]
    shape=tuple(d+1 for d in degrees); A=np.empty(shape,dtype=object);A.fill(Q(0))
    for e,v in D.items(): A[e]=v
    for axis,(lo,hi) in enumerate(box):
        d=degrees[axis];n=d+1;delta=hi-lo
        BT=Mat([[sum((Q(math.comb(e,j)*math.comb(i,j),math.comb(d,j))*delta**j*lo**(e-j) for j in range(min(i,e)+1)),Q(0)) for e in range(n)] for i in range(n)])
        V=np.moveaxis(A,axis,0); sh=V.shape; cols=math.prod(sh[1:]); M=BT*Mat(V.reshape(n,cols).tolist())
        A=np.moveaxis(np.array([[M[i,j] for j in range(cols)] for i in range(n)],dtype=object).reshape(sh),0,axis)
    lo=np.empty(shape);hi=np.empty(shape)
    for e in np.ndindex(shape):
        z=arb(A[e]);lo[e]=np.nextafter(float(z.lower()),-np.inf);hi[e]=np.nextafter(float(z.upper()),np.inf)
    assert np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))
    return lo,hi


def verify_box(D,box,far=False,limit=500000):
    lo,hi=bernstein(D,box)
    stack=[(lo,hi,((0.,1.),)*3,0)]; checked=passed=excluded=maxdepth=0
    while stack:
        l,h,unit,depth=stack.pop();checked+=1;maxdepth=max(maxdepth,depth)
        if np.max(h)<0: passed+=1;continue
        physical=[(arb(a)+arb(b-a)*arb(c),arb(a)+arb(b-a)*arb(d)) for (a,b),(c,d) in zip(box,unit)]
        if far and physical[2][0]>(physical[0][1].sqrt()+physical[1][1].sqrt())**2:
            excluded+=1;continue
        # A strictly positive value at any included corner is a counterexample.
        for idx in itertools.product((0,-1),repeat=3):
            if l[idx]>0:
                coords=[physical[k][0 if idx[k]==0 else 1] for k in range(3)]
                if not far or coords[2]<=(coords[0].sqrt()+coords[1].sqrt())**2:
                    raise AssertionError(('POSITIVE_CORNER',unit,idx,float(l[idx]),list(map(str,coords))))
        if checked>=limit or depth>=65:
            raise AssertionError(('BOX_LIMIT',unit,depth,float(np.max(h)),len(stack)))
        mid=(l+h)*.5
        scores=[float(np.max(np.abs(np.diff(mid,axis=i)))) if mid.shape[i]>1 else -1 for i in range(3)]
        axis=int(np.argmax(scores));ll,lh,rl,rh=split(l,h,axis)
        a,b=unit[axis];m=(a+b)*.5;bx1=list(unit);bx2=list(unit);bx1[axis]=(a,m);bx2[axis]=(m,b)
        stack.extend([(rl,rh,tuple(bx2),depth+1),(ll,lh,tuple(bx1),depth+1)])
    return {'checked':checked,'passed':passed,'excluded':excluded,'depth':maxdepth}


def verify(path,only=None):
    t0=time.time();r=json.load(open(path));n=int(r['dimension']);K=int(r['radial_nodes']);L=int(r['angular_degree']);D=int(r['radial_degree'])
    r2=Q(str(Fraction(r['exclusion_squared'])));R2=Q(str(Fraction(r['truncation_squared'])));assert 3<=n<=48 and 0<r2<R2<2*r2
    As=[[[Q(z) for z in row] for row in A] for A in r['matrices']]
    assert len(As)==L+1
    for k,A in enumerate(As):
        m=K+1 if k==0 else K;assert len(A)==m and all(len(row)==m for row in A)
        assert all(A[i][j]==A[j][i] for i in range(m) for j in range(m))
        B=[row[:] for row in A]
        for j in range(m):
            assert B[j][j]>0,('PSD',k,j,str(B[j][j]))
            for i in range(j+1,m):
                for h in range(j+1,m): B[i][h]-=B[i][j]*B[j][h]/B[j][j]
    print('Exact rational PSD: PASS',flush=True)
    X=Poly([0,1]);P=Poly(list(map(Q,r['physical_polynomial'])));H=Poly(list(map(Q,r['fourier_polynomial'])));a=Q(n,2)
    ls=[Poly([1]),a-2*X]
    for k in range(1,D):ls.append(((2*k+a-2*X)*ls[-1]-(k+a-1)*ls[-2])/Q(k+1))
    rem=P;fh=Poly([])
    for k in range(D,-1,-1):
        c=rem[k]/ls[k][k];rem-=c*ls[k];fh+=(-1)**k*c*ls[k]
    assert not rem and fh==H and H[0]==Q(r['normalization'])
    def variation(vals):
        ss=[1 if v>0 else -1 for v in vals if v];return sum(a!=b for a,b in zip(ss,ss[1:]))
    def roots(p):
        seq=[p,p.derivative()]
        while seq[-1].degree()>0:
            q=-(seq[-2]%seq[-1])
            if not q:break
            num=q.numer();seq.append(Poly(num/abs(num.content())))
        return variation([q[0] for q in seq])-variation([q[q.degree()] for q in seq])
    pi_lo=Q('314159265358979323846264338327950288419716939937510/100000000000000000000000000000000000000000000000000')
    assert arb(pi_lo)<arb.pi()
    assert H[0]>0
    physical_origin_correction,corrected_norm=verify_guards(r,P,H,r2,R2,pi_lo)
    print('Exact Fourier identity and reciprocal-guard global signs: PASS',flush=True)
    intervals=[(r2+(R2-r2)*i/(K-1),r2+(R2-r2)*(i+1)/(K-1)) for i in range(K-1)]
    hats=[]
    for i,(lo,hi) in enumerate(intervals):
        w=[Poly([]) for _ in range(K)];w[i]=(hi-X)/(hi-lo);w[i+1]=(X-lo)/(hi-lo);hats.append(w)
    if only in (None,'pair'):
        P1=P.derivative()-P;P2=P.derivative().derivative()-2*P.derivative()+P;pi=arb.pi();checked=0;maxdepth=0
        for cell,(lo,hi) in enumerate(intervals):
            w=hats[cell];q=sum((2*As[0][0][j+1]*w[j] for j in range(K)),Poly([]))
            for ell,A in enumerate(As):
                off=int(ell==0)
                radial=sum((A[j+off][k+off]*w[j]*w[k] for j in range(K) for k in range(K)),Poly([]))
                q+=radial*(X/r2)**ell
            stack=[(lo,hi,0)]
            while stack:
                lo,hi,d=stack.pop();mid=(lo+hi)/2;hh=(hi-lo)/2;u=arb(mid);h=arb(hh);xx=pi*u;ex=(-xx).exp()
                g=ev(P,xx)*ex+ev(q,u);gp=pi*ev(P1,xx)*ex+ev(q.derivative(),u)
                M2=pi*pi*mag(P2,xx,pi*h)*(-pi*arb(lo)).exp()+mag(q.derivative().derivative(),u,h)
                upper=g+abs(gp)*h+M2*h*h/2;checked+=1;maxdepth=max(maxdepth,d)
                if upper<0:continue
                assert not g>0,('POSITIVE_PAIR_MIDPOINT',str(mid),str(g))
                assert d<40 and checked<100000,('PAIR_LIMIT',str(lo),str(hi),str(upper))
                stack.extend([(lo,mid,d+1),(mid,hi,d+1)])
        print('Global solid-harmonic pair inequality: PASS',checked,maxdepth,flush=True)
    one={(0,0,0):Q(1)};lin={(1,0,0):1/(2*r2),(0,1,0):1/(2*r2),(0,0,1):-1/(2*r2)}
    polys=[one,lin]
    for ell in range(1,L):
        polys.append(add(mul(lin,polys[-1]),{},Q(1)))
        polys[-1]={e:v*Q(2*ell+n-2,ell+n-2) for e,v in polys[-1].items()}
        polys[-1]=add(polys[-1],shift(polys[-3],(1,1,0)),-Q(ell,ell+n-2)/r2**2)
    cache={}
    def kernel(i,j):
        if (i,j) in cache:return cache[i,j]
        out={}
        for ell,A in enumerate(As):
            off=int(ell==0);rad={}
            for b in (0,1):
                for c in (0,1):
                    z=sum((hats[i][ii][b]*A[ii+off][jj+off]*hats[j][jj][c] for ii in range(K) for jj in range(K)),Q(0))
                    if z:rad[b,c,0]=z
            out=add(out,mul(rad,polys[ell]))
        cache[i,j]=out;return out
    report=[]
    if only in (None,'cyclic'):
        for i,j,k in itertools.combinations_with_replacement(range(K-1),3):
            poly=add(add(kernel(i,j),permute(kernel(i,k),(0,2,1))),permute(kernel(j,k),(2,0,1)))
            out=verify_box(poly,[intervals[i],intervals[j],intervals[k]])
            report.append(('cyclic',i,j,k,out));print('CYCLIC',i,j,k,out,'seconds',time.time()-t0,flush=True)
        print('Entire compact cyclic domain: PASS',flush=True)
    if only in (None,'far'):
        for i,j in itertools.combinations_with_replacement(range(K-1),2):
            out=verify_box(kernel(i,j),[intervals[i],intervals[j],(R2,4*R2)],far=True)
            report.append(('far',i,j,out));print('FAR',i,j,out,'seconds',time.time()-t0,flush=True)
        print('Entire far domain: PASS',flush=True)
    objective=(P[0]+physical_origin_correction+As[0][0][0])/corrected_norm;assert objective==Q(r['exact_objective'])
    density=arb.pi()**arb(a)/(arb(a)+1).gamma()*arb(r2/4)**arb(a)*arb(objective)
    print('Density enclosure:',density,flush=True)
    if only is None:
        assert density<arb(r['certified_density_upper'])
        print('PASS: complete reciprocal-Gaussian-guarded cyclic certificate.',flush=True)
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('path');ap.add_argument('--only',choices=['pair','cyclic','far']);args=ap.parse_args();verify(args.path,args.only)
