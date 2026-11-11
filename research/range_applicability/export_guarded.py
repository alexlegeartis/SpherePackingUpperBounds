"""Discovery/export only; every generated file must pass verify_guarded.py."""
import argparse,json,math
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
from reciprocal_guards import bernstein_upper,certify_complement
ctx.prec=256
PI_LO=Q('314159265358979323846264338327950288419716939937510/100000000000000000000000000000000000000000000000000')


def islands(p,start,desired_sign):
    real=[]
    for z,m in p.complex_roots():
        if z.imag==0 and float(z.real)>float(start):real.append(float(z.real))
    real.sort();points=[float(start)]+real;bad=[]
    for i in range(len(points)-1):
        mid=Q(str(Fraction(str((points[i]+points[i+1])/2))))
        if desired_sign*p(mid)<0:
            lo=start if i==0 else Q(math.floor(points[i]*1000)-1,1000)
            hi=Q(math.ceil(points[i+1]*1000)+1,1000)
            bad.append((max(start,lo),hi))
    merged=[]
    for lo,hi in bad:
        if merged and lo<=merged[-1][1]:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
        else:merged.append((lo,hi))
    certify_complement(p,start,merged,desired_sign)
    return merged


def outward(x):
    if x<=0:return Q(0)
    value=float(x.upper());exp=math.floor(math.log10(value))-1
    step=Q(10)**exp;q=Q(math.ceil(value/float(step)))*step
    while not arb(q)>x:q+=step
    return q


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('output');ap.add_argument('--upper',required=True);args=ap.parse_args()
    r=json.load(open(args.source));n=r['dimension'];s=Q(str(Fraction(r['exclusion_squared'])));T=Q(str(Fraction(r['truncation_squared'])));P=Poly(list(map(Q,r['physical_polynomial'])));H=Poly(list(map(Q,r['fourier_polynomial'])))
    c=Q(math.floor(float(arb.pi()*arb(s))*100),100);assert arb(c)<arb.pi()*arb(s)
    phys=islands(P,T*PI_LO,-1);y0=(2*n-c)/16;freq=islands(H,y0,1)
    eta=Q(0);zeta=Q(0)
    for lo,hi in phys:
        req=arb(bernstein_upper(P,lo,hi))*(-arb(3)*arb(lo)/4).exp()/arb(lo-c)
        eta=max(eta,outward(req))
    for lo,hi in freq:
        req=arb(2**n)*arb(bernstein_upper(-H,lo,hi))*(-arb(3)*arb(lo)/4).exp()/arb(c-Q(n,8)+lo/16)
        zeta=max(zeta,outward(req))
    D=Q(2**n)*(2*n-c);E=(c-Q(n,8))/Q(2**n);norm=H[0]-eta*D+zeta*E
    assert norm>0
    obj=(P[0]+c*(eta+zeta)+Q(r['matrices'][0][0][0]))/norm
    r['base_exact_objective']=r['exact_objective'];r['exact_objective']=str(obj);r['kernel_family']='solid_harmonic_cyclic_reciprocal_guards';r['reciprocal_gaussian_guards']={'c':str(c),'eta':str(eta),'zeta':str(zeta),'physical_intervals':[[str(a),str(b)]for a,b in phys],'fourier_intervals':[[str(a),str(b)]for a,b in freq],'corrected_normalization':str(norm)};r['certified_density_upper']=args.upper;r['certification_status']='UNVERIFIED reciprocal Gaussian guards; requires complete verify_guarded.py replay.'
    Path(args.output).write_text(json.dumps(r,indent=2)+'\n')
    a=arb(n)/2;density=arb.pi()**a/(a+1).gamma()*arb(s/4)**a*arb(obj)
    print('Guard parameters',r['reciprocal_gaussian_guards'],flush=True);print('Unverified guarded objective density',density,flush=True)
