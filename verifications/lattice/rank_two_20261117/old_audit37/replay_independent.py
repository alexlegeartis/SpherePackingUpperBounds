if not __debug__:raise RuntimeError('Assertions required')
import sys,json,math,hashlib
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
from flint import arb,ctx
sys.set_int_max_str_digits(0);ctx.prec=1024
src=Path(sys.argv[1]);raw=src.read_bytes();D=json.loads(raw);n=D['n'];p=list(map(F,D['p']));q=list(map(F,D['q']));ds=list(map(F,D['d']));K=len(p)-1;tau=F(D['tau']);E=list(map(F,D['E']));g=list(map(F,D['g']));J=E.index(F(2));M=len(g);nn=M-J;C=F(D['C']);origin_ratio=F(D['origin_ratio']);assert F(D['origin_C'])==origin_ratio**2*C and origin_ratio>=0;angular=D['cyclic_angular'];B=F(angular['B_exact']);qm=sum(F(v)*(-1)**i for i,v in enumerate(angular['gegenbauer_coefficients']));hpow=-int(angular['profile_a_exact']);beta=F(angular['profile_b_exact']);assert n==37 and C>=0 and E[0]==1 and E[-1]==F(9,4) and len(E)==M+1 and min(g)>=0
assert angular==json.load(open('cyclic/angular.json'))
for j in range(K+1):
 den=math.factorial(j)*math.prod(F(n,2)+a for a in range(j));values=[ds[k]*F(math.factorial(k),math.factorial(k-j)) for k in range(j,K+1)]
 assert p[j]==(-1)**j*sum(values)/den
 assert q[j]==(-1)**j*sum(v*(-1)**k for k,v in zip(range(j,K+1),values))/den
print('polynomial reconstruction PASS',flush=True)
# independent integer Bernstein global replay
for name,poly,start in [('q',q,F(0)),('minus_p',[-v for v in p],tau*E[-1])]:
 vals=[sum(poly[k]*math.comb(k,j)*start**(k-j) for k in range(j,K+1))/math.comb(K,j) for j in range(K+1)];den=math.lcm(*(x.denominator for x in vals));pending=[([int(v*den) for v in vals],0)];leaves=0
 while pending:
  b,depth=pending.pop()
  if min(b)>0:leaves+=1;continue
  assert b[0]>0 and b[-1]>0 and depth<100
  left=[b[0]<<K];right=[b[-1]<<K]
  for level in range(1,K+1):b=[u+v for u,v in zip(b,b[1:])];left.append(b[0]<<(K-level));right.append(b[-1]<<(K-level))
  for row in (left,right[::-1]):
   div=math.gcd(*row);pending.append(([v//div for v in row],depth+1))
 print(name,'global PASS',leaves,flush=True)
def ar(v):return arb(v.numerator)/v.denominator
def translated(coeff,a):return [sum(coeff[r]*math.comb(r,j)*a**(r-j) for r in range(j,len(coeff))) for j in range(len(coeff))]
def radial_bound(left,right):
 order=18;width=right-left;Fpoly=translated(p,tau*left);Fpoly=[v*tau**j for j,v in enumerate(Fpoly)]
 a=1+qm/B
 terms=[(F(1),tau/2,F(0),Fpoly,None),(C*a/B,2*beta,2*beta,None,2*hpow),(-2*C*origin_ratio/B,beta,beta,None,hpow),(2*C*4**hpow/B,5*beta,2*beta,None,2*hpow)]
 total=[arb(0) for _ in range(order+1)];remainder=arb(0)
 for multiplier,rate,constant,pol,degree in terms:
  if pol is None:pol=[F(math.comb(degree,j))*left**(degree-j) for j in range(degree+1)]
  amplitude=ar(multiplier)*(ar(constant-rate*left)).exp();absol=ar(abs(multiplier))*(ar(constant-rate*left)).exp()
  for k in range(order+1):total[k]+=amplitude*ar(sum(pol[j]*(-rate)**(k-j)/math.factorial(k-j) for j in range(min(k,len(pol)-1)+1)))
  k=order+1;bound=F(0)
  for j in range(min(k,len(pol)-1)+1):
   sup=(F(math.comb(degree,j))*right**(degree-j) if degree is not None else sum(abs(pol[r])*math.comb(r,j)*width**(r-j) for r in range(j,len(pol))))
   bound+=sup*rate**(k-j)/math.factorial(k-j)
  remainder+=absol*ar(bound*width**k)
 return [sum((total[j]*ar(width**j*F(math.comb(i,j),math.comb(order,j))) for j in range(i+1)),arb(0))+remainder for i in range(order+1)]
for i,paths in enumerate(D['radial_leaves']):
 cells=[]
 for path in paths:
  assert set(path)<=set('01');den=2**len(path);num=int(path,2) if path else 0;cells.append((F(num,den),F(num+1,den)))
 cells.sort();assert cells[0][0]==0 and cells[-1][1]==1 and all(a[1]==b[0] for a,b in zip(cells,cells[1:]))
 for lo,hi in cells:
  left=E[i]+(E[i+1]-E[i])*lo;right=E[i]+(E[i+1]-E[i])*hi
  assert all(v<=ar(g[i]) for v in radial_bound(left,right)),('radial',i,lo,hi)
print('all mixed radial intervals PASS',M,flush=True)
report=dict(status='PASS',sha256=hashlib.sha256(raw).hexdigest(),global_signs=True,radial_intervals=M,radial_cells=sum(map(len,D['radial_leaves'])),independent_taylor_order=18,precision=1024)
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
