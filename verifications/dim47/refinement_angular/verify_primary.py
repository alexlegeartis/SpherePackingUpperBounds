"""Complete frozen 47D cyclic density certificate replay. No optimizer."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,subprocess,hashlib
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Do not use Python -O')
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).resolve().parent;S=D/'continuum_three_point'
expected={'cyclic_angular40_exact.json': 'c2e67ffcd6e5c94a33d94d9276108fd5955b45c5', 'independent_cyclic40_angular.py': '8f73359c650845b833c5372d47014bda737cd51c', 'independent_cyclic40_domain.py': '2dceea140bfd317757c68e46792cfa6f50f2925a'}
for name,want in expected.items():
 b=(S/name).read_bytes();assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==want
for part in ['angular','domain']:
 result=subprocess.run([sys.executable,str(S/f'independent_cyclic40_{part}.py')],cwd=D,capture_output=True,text=True)
 (D/f'source_{part}_replay.log').write_text(result.stdout+result.stderr);print(result.stdout,flush=True);assert result.returncode==0
r=json.loads((D/'certificate.json').read_text());ang=json.loads((S/'cyclic_angular40_exact.json').read_text());assert r['dimension']==47
p=list(map(F,r['physical_polynomial']));h=list(map(F,r['fourier_polynomial']));M=F(r['M_exact']);B=F(ang['B_exact']);n=47;a=F(47,2)
assert len(p)==len(h)==112 and h[0]>0 and M>0
rf=lambda a,k:math.prod((a+j for j in range(k)),start=F(1))
q=lambda x:fmpq(x.numerator,x.denominator);af=lambda x:arb(q(x));ff=lambda x:F(int(x.numerator),int(x.denominator))
transform=[sum((p[k]*math.factorial(k)*(-1)**j*rf(a+j,k-j)/math.factorial(k-j)/math.factorial(j) for k in range(j,len(p))),F(0)) for j in range(len(p))]
assert transform==h
Q=fmpq_poly(list(map(fmpq,ang['polynomial_coefficients'])));t=fmpq_poly([0,1]);G=[fmpq_poly([1]),t]
for j in range(1,Q.degree()):G.append(((2*j+n-2)*t*G[-1]-j*G[-2])/(j+n-2))
rem=Q;cs=[F(0)]*len(G)
for j in range(len(G)-1,-1,-1):cs[j]=ff(rem[j]/G[j][j]);rem-=q(cs[j])*G[j]
assert not rem and min(cs)>=0 and all(g(1)==1 for g in G) and sum(cs)==B
b0=cs[0];assert b0>0 and b0==F(r['b0_exact'])
assert b0==sum((ff(Q[2*j])*rf(F(1,2),j)/rf(a,j) for j in range(Q.degree()//2+1)),F(0))
origin=M*B/b0;degree0=M*b0/B;assert origin>0 and degree0>0 and origin*degree0==M*M
P=fmpq_poly(list(map(q,p)));H=fmpq_poly(list(map(q,h)))
def ab(x):
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(110));return lo,lo+x**221/221
al,au=ab(F(1,5));bl,bu=ab(F(1,239));pl=16*al-4*bu;pu=16*au-4*bl;assert F(51,4)*pl>40
for name,poly,cut in [('fourier',H,F(0)),('physical',P,7*pl)]:
 roots=poly.numer().complex_roots();assert sum(m for z,m in roots)==111
 assert all(not z.imag.contains(0) or z.real<af(cut) for z,m in roots)
 print(name,'global sign root isolation PASS',flush=True)
assert H(0)>0 and P(q(7*pl))<0
AP=arb_poly([af(v) for v in p]);stack=[(6*pl,7*pu,0)];passed=0;maxdepth=0
while stack:
 l,u,dep=stack.pop();mid=(l+u)/2;rad=(u-l)/2
 co=AP(arb_poly([af(mid),af(rad)]));up=co[0]+sum((abs(co[j]) for j in range(1,len(co))),arb(0));ok=up<0
 if not ok:
  cap=2-(-arb(fmpq(17,8))*(af(l)-6*af(pu))).exp()
  expupper=(arb(fmpq(9,8))*af(u)-arb(fmpq(51,4))*af(pl)).exp()
  ok=up*expupper<af(M)*cap
 if ok:passed+=1;maxdepth=max(maxdepth,dep);continue
 assert dep<32
 stack.extend([(mid,u,dep+1),(l,mid,dep+1)])
obj=(p[0]+origin)/h[0];assert obj==F(r['objective_exact']);upper=F('382349047/625000000000000')
rational_upper_squared=obj**2*4*pu**46*6**47*math.factorial(24)**2/math.factorial(48)**2;assert obj>0 and 0<rational_upper_squared<upper**2
baseline=F('0.0000006490757338');assert upper<baseline
density=af(obj)*(arb.pi()*arb(3)/2)**(arb(47)/2)/(arb(49)/2).gamma();assert density<af(upper)
out={'pass':True,'dimension':47,'strict_upper':'0.0000006117584752','density_interval':str(density),'radial_Fourier_identity':True,'global_Fourier_and_physical_signs':True,'coupling_intervals':passed,'coupling_depth':maxdepth,'angular_dimension_specific_PSD':True,'origin_rank_one_Gram':True,'new_complete_cyclic_domain':True,'exact_rational_rounding':True,'pinned_comparator':'0.0000006490757338','fractional_improvement_lower':str(1-upper/baseline)}
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
