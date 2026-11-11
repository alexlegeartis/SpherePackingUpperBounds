"""Independent27D rational Fourier, Sturm, harmonic PSD and coupled proof."""
import json,math,sys
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq, fmpq_poly, arb,arb_poly,ctx
assert __debug__
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).resolve().parent;d=json.loads((D/'certificate.json').read_text());angular=json.loads((D/'angular.json').read_text());assert d['dimension']==27
assert F(d['r2'])==F(27,8) and F(d['R2'])==7 and F(d['kappa'])==F(34,9)
p=list(map(F,d['physical_polynomial']));h=list(map(F,d['fourier_polynomial']));assert len(p)==len(h)==112
M=F(d['M_exact']);N=h[0];assert M>0 and N>0
q=lambda x:fmpq(x.numerator,x.denominator)
a=lambda x:arb(q(x))
def rise(z,k):return math.prod((z+i for i in range(k)),start=F(1))
# Direct monomial Fourier transform in27D (Laguerre parameter25/2).
fh=[sum((p[k]*math.factorial(k)*(-1)**j*rise(F(27,2)+j,k-j)/math.factorial(k-j)/math.factorial(j) for k in range(j,112)),F(0))for j in range(112)]
assert fh==h
pol=list(map(F,angular['polynomial_coefficients']));assert len(pol)==20 and not any(pol[16:]);B=F(angular['B_exact']);assert len(pol)==20 and sum(pol)==B
basis=[]
for ell in range(20):
 v=[F(0)]*20;norm=rise(F(25),ell)/math.factorial(ell)
 for j in range(ell//2+1):v[ell-2*j]=(-1)**j*rise(F(25,2),ell-j)*2**(ell-2*j)/(math.factorial(j)*math.factorial(ell-2*j)*norm)
 assert sum(v)==1;basis.append(v)
rem=pol[:];cs=[F(0)]*20
for ell in range(19,-1,-1):
 cs[ell]=rem[ell]/basis[ell][ell]
 for j in range(20):rem[j]-=cs[ell]*basis[ell][j]
assert not any(rem) and min(cs)>=0 and sum(cs)==B
b0=cs[0];assert b0==F(d['b0_exact']) and b0==sum((pol[2*j]*rise(F(1,2),j)/rise(F(27,2),j)for j in range(10)),F(0))
origin=M*B/b0;diag=M*b0/B;assert origin>0 and diag>0 and origin*diag==M*M

def atan_bounds(x):
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1)for j in range(100));return lo,lo+x**201/201
l5,u5=atan_bounds(F(1,5));l239,u239=atan_bounds(F(1,239));pl=16*l5-4*u239;pu=16*u5-4*l239
assert F(51,4)*pl>40
P=fmpq_poly(list(map(q,p)));H=fmpq_poly(list(map(q,h)))
def sturm(poly):
 seq=[poly,poly.derivative()]
 while seq[-1].degree()>0:
  r=-(seq[-2]%seq[-1]);assert r
  z=r.numer();seq.append(fmpq_poly(z/abs(z.content())))
 return seq
def variation(seq,t):
 vals=[v[v.degree()] if t is None else v(q(t))for v in seq];sg=[int(v>0)-int(v<0)for v in vals if v]
 return sum(x!=y for x,y in zip(sg,sg[1:]))
for name,poly,cut,sign in [('Fourier',H,F(0),1),('physical',P,7*pl,-1)]:
 print('Independent exact Sturm start',name,flush=True);s=sturm(poly);v0=variation(s,cut);vi=variation(s,None)
 assert v0==vi and sign*poly(q(cut))>0;print('Independent exact Sturm PASS',name,v0,vi,flush=True)
# Independent second-derivative enclosure of the full divided coupling function.
AP=arb_poly([a(x)for x in p]);beta=arb(25)/9;k=arb(34)/9;C=51*arb.pi()/4;Am=a(M)
DP=AP.derivative()+beta*AP;D2P=AP.derivative().derivative()+2*beta*AP.derivative()+beta**2*AP
stack=[(F(27,8)*pl,7*pu,0)];leaves=0;depthmax=0
while stack:
 lo,hi,depth=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2;v=a(mid);rr=a(rad);ea=(beta*v-C).exp();eb=(-k*v+C).exp()
 center=AP(v)*ea-2*Am+Am*eb;first=DP(v)*ea-k*Am*eb
 taylor=D2P(arb_poly([v,rr]));mag=sum((abs(x)for x in taylor.coeffs()),arb(0))
 second=mag*(beta*a(hi)-C).exp()+k*k*Am*(-k*a(lo)+C).exp()
 upper=center+abs(first)*rr+second*rr*rr/2
 if upper<0:leaves+=1;depthmax=max(depthmax,depth);continue
 assert depth<35;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
objective=(p[0]+origin)/N;assert objective==F(d['objective_exact']) and objective>0
threshold=F('0.0006864321082');rational_upper=objective**2*4*pu**26*F(27,8)**27*math.factorial(14)**2/math.factorial(28)**2;assert rational_upper<threshold**2
volume=2*arb.pi()**13*a(F(27,8))**13*a(F(27,8)).sqrt()*math.factorial(14)/math.factorial(28);density=a(objective)*volume
out={'pass':True,'dimension':27,'Fourier_identity':True,'global_exact_Sturm_signs':True,'Gegenbauer_PSD':True,'origin_Gram_PSD':True,'pair_leaves':leaves,'pair_depth':depthmax,'density_interval':str(density),'strict_upper':str(threshold),'strict_upper_decimal':'0.0006864321082','inherited_dependency':'New27D angular and global cyclic verification remains required'}
(D/'audit.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
