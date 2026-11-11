import sys,json,math
sys.set_int_max_str_digits(0);sys.path.insert(0,'solver_packages')
from fractions import Fraction as F
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
from pathlib import Path
ctx.prec=512;p=Path('continuum_three_point');r=json.loads((p/'adaptive_cyclic40_radial111.json').read_text());ang=json.loads((p/'cyclic_angular40_exact.json').read_text());B=F(ang['B_exact']);assert B==F(r['B']);assert F(str(r['kappa']))==F(17,8)
ps=list(map(F,r['physical_polynomial']));hs=list(map(F,r['fourier_polynomial']));norm=hs[0];assert norm>0;ps=[x/norm for x in ps];hs=[x/norm for x in hs]
qs=[sum(ps[k]*math.factorial(k)*(-1)**j*F(math.comb(k+23,k-j),math.factorial(j)) for k in range(j,len(ps))) for j in range(len(ps))];assert qs==hs and hs[0]==1
ff=lambda x:fmpq(x.numerator,x.denominator);P=fmpq_poly(list(map(ff,ps)));H=fmpq_poly(list(map(ff,hs)));M=F(r['neighbor_cost'])/(B*norm)*F(1000001,1000000);assert M>0
# Independent tighter Machin bounds, 110 alternating terms.
def ab(x):
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(110));return lo,lo+x**221/221
al,au=ab(F(1,5));bl,bu=ab(F(1,239));pl=16*al-4*bu;pu=16*au-4*bl;assert F(51,4)*pl>40
for poly,cut in [(H,F(0)),(P,7*pl)]:
 roots=poly.numer().complex_roots();assert sum(m for z,m in roots)==111;assert all(not z.imag.contains(0) or z.real<arb(ff(cut)) for z,m in roots)
assert P(ff(7*pl))<0
AP=arb_poly([arb(ff(q)) for q in ps]);stack=[(6*pl,7*pu,0)];passed=0;maxdepth=0
while stack:
 a,b,depth=stack.pop();mid=(a+b)/2;rad=(b-a)/2
 coeff=AP(arb_poly([arb(ff(mid)),arb(ff(rad))]));up=coeff[0]+sum((abs(coeff[i]) for i in range(1,len(coeff))),arb(0))
 # On the physical domain s>=6*pi, division by h leaves RHS>=M>0.
 ok=up<0
 if not ok:
  cap=2-(-arb(fmpq(17,8))*(arb(ff(a))-6*arb(ff(pu)))).exp()
  expupper=(arb(fmpq(9,8))*arb(ff(b))-arb(fmpq(51,4))*arb(ff(pl))).exp()
  ok=up*expupper<arb(ff(M))*cap
 if ok:passed+=1;maxdepth=max(maxdepth,depth);continue
 assert depth<30
 stack.extend([(mid,b,depth+1),(a,mid,depth+1)])
U=ps[0]+B*M;assert U<F(17835,1000);ref=json.loads((p/'cyclic40_radial111_exact.json').read_text());assert U==F(ref['objective_exact']) and M==F(ref['M_exact'])
density=F(17835,1000)*(F(355,113)*F(3,2))**24/math.factorial(24)
out={'pass':True,'all112_Fourier_coefficients_match':True,'precision':512,'Machin_terms':110,'cap_intervals':passed,'max_depth':maxdepth,'objective_exact':str(U),'objective_float':float(U),'strict_upper':'17835/1000','density_upper':str(density),'density_upper_float':float(density)}
(p/'independent_cyclic40_radial.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k not in ('objective_exact','density_upper')})
