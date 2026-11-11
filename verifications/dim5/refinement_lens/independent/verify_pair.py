"""Independent full6D pair inequality, with exact count-threshold boundaries."""
import json
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
assert __debug__
ctx.prec=384;D=Path(__file__).resolve().parent;s=json.loads((D/'witness.json').read_text());raw=json.loads((D/'witness_psd_exact.json').read_text());rad=json.loads((D/'radial_exact.json').read_text())
q=lambda x:Q(x.numerator,x.denominator)
a=lambda x:arb(q(x))
As=[[[Q(z)for z in row]for row in A]for A in raw['matrices']];R=F(s['R']);X=Poly([0,1]);z=(2*X-q(R+1))/q(R-1);T=[Poly([1]),z]
for j in range(1,4):T.append(2*z*T[-1]-T[-2])
Kpair=2*sum((As[0][0][i+1]*T[i]for i in range(5)),Poly([]))
for ell,A in enumerate(As):
 off=int(ell==0)
 Kpair+=sum((A[i+off][j+off]*T[i]*T[j]for i in range(5)for j in range(5)),Poly([]))
AP=arb_poly([arb(Q(x))for x in rad['P']]);KP=arb_poly([arb(x)for x in Kpair.coeffs()]);P1=AP.derivative()-AP;P2=AP.derivative().derivative()-2*AP.derivative()+AP
radii=list(map(F,s['count_radii']));weights=[F(str(x))for x in s['count_weights']];edges=sorted(set([F(1),R]+radii));stack=[(lo,hi,0)for lo,hi in zip(edges[:-1],edges[1:])];leaves=maxdepth=0;pi=arb.pi()
while stack:
 lo,hi,depth=stack.pop();mid=(lo+hi)/2;radius=(hi-lo)/2;r=a(mid);rr=a(radius);xp=pi*r*r
 allowance=sum((w for threshold,w in zip(radii,weights)if mid<threshold),F(0))
 value=AP(xp)*(-xp).exp()+KP(r)-a(allowance)
 derivative=2*pi*r*P1(xp)*(-xp).exp()+KP.derivative()(r)
 xm=pi*a((lo*lo+hi*hi)/2);xr=pi*a((hi*hi-lo*lo)/2)
 mag1=sum((abs(x)for x in P1(arb_poly([xm,xr])).coeffs()),arb(0));mag2=sum((abs(x)for x in P2(arb_poly([xm,xr])).coeffs()),arb(0))
 second=(2*pi*mag1+4*pi*pi*a(hi*hi)*mag2)*(-pi*a(lo*lo)).exp()+abs(KP.derivative().derivative()(a(lo).union(a(hi))))
 upper=value+abs(derivative)*rr+second*rr*rr/2
 if upper<0:leaves+=1;maxdepth=max(depth,maxdepth);continue
 assert depth<30;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
# At each threshold, the strict count indicator removes that threshold's weight.
for threshold in radii:
 r=a(threshold);xp=pi*r*r;allowance=sum((w for t,w in zip(radii,weights)if threshold<t),F(0))
 assert AP(xp)*(-xp).exp()+KP(r)<a(allowance)
out={'pass':True,'closed_interval_leaves':leaves,'maxdepth':maxdepth,'thresholds_checked_with_smaller_allowance':[str(r)for r in radii],'scope':'Base pair passes; Gaussian guard is nonpositive for every r>=1'}
(D/'pair_audit.json').write_text(json.dumps(out,indent=2));print(out)
