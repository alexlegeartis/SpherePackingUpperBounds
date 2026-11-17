import json,math,sys,hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
ctx.prec=1024;sys.set_int_max_str_digits(0)
pth=Path(__file__).parent;old=json.loads((pth.parent/'certificate36.json').read_text());cand=json.loads((pth/'input_candidate.json').read_text());basis=json.loads((pth/'input_basis.json').read_text());tau=F(basis['tau']);qzero=F(basis['qzero']);qnew=F(basis['q']);polys=[list(map(F,v))for v in basis['polys']];assert tau==F(old['tau']) and qzero==F(old['q'][0]);assert polys[0]==list(map(F,old['p']))
import sympy as sp
x=sp.symbols('x');co=[sp.Rational(str(v))for v in json.loads((pth.parent/'rank2_candidate37_B17677213_20261117.json').read_text())['even_coefficients']];Q=sp.Poly(sum(v*sp.gegenbauer(2*j,sp.Rational(35,2),x)/sp.gegenbauer(2*j,sp.Rational(35,2),1)for j,v in enumerate(co))-sp.Rational(2,10000),x);rest=Q;target=[]
for j in range(Q.degree()//2,-1,-1):
 base=sp.Poly(sp.gegenbauer(2*j,17,x)/sp.gegenbauer(2*j,17,1),x);c=rest.nth(2*j)/base.nth(2*j);assert c>=0;target.append(c);rest-=c*base
assert rest.is_zero and target[-1]>0 and qnew==F(Q.eval(1)/target[-1])

def op(p):return [((36+4*j)*(p[j]if j<len(p)else 0)-(p[j-1]if j>=1 else 0)-(j+1)*(72+4*j)*(p[j+1]if j+1<len(p)else 0))/tau for j in range(len(p)+1)]
assert len(polys)==3
assert op(polys[0])==polys[1] and op(polys[1])==polys[2]
v=F(str(cand['v']));w=max(F(str(cand['w'])),v*v/4+F(1,10**16));assert w>0 and w-v*v/4>0
p=[sum(z*(poly[j]if j<len(poly)else 0)for z,poly in zip([F(1),v,w],polys))for j in range(max(map(len,polys)))];K=len(p)-1
Z=[[F(str(z))for z in row]for row in cand['origin']];assert Z[0][1]==Z[1][0];z=Z[0][1]-F(1,10**10);aZ=Z[1][1];assert aZ>0 and z<0;s=z*z/aZ+F(1,10**20);assert s*aZ-z*z>0
rho=F(1,10**70);assert F(168**168,math.factorial(168))>1/rho;factor=(1+rho)/(1-rho);A=2*qzero*aZ*factor/qnew;B=2*qzero*z/qnew;assert A>0 and B<0 and A+B<0
E=list(map(F,old['E']));g=[max(F(0),F(str(x)))+F(1,10**24)for x in cand['g']];M=len(g);J=E.index(F(2));nn=M-J;assert len(E)==M+1
ds=list(map(F,old['d']));oldp=list(map(F,old['p']));oldq=list(map(F,old['q']));N=len(ds)-1
for j in range(N+1):
 den=math.factorial(j)*math.prod(F(18)+a for a in range(j));values=[ds[k]*F(math.factorial(k),math.factorial(k-j))for k in range(j,N+1)]
 assert oldp[j]==(-1)**j*sum(values)/den
 assert oldq[j]==(-1)**j*sum(a*(-1)**k for k,a in zip(range(j,N+1),values))/den
assert hashlib.sha256((pth.parent/'rank2_candidate37_B17677213_20261117.json').read_bytes()).hexdigest()=='a57da6cbe33f3f7395ef40b76d3450388115423d60f83270920b10fc4f97212d'
print('exact multiplier basis, PSD blocks, radial tail coupling and rho PASS',flush=True)
def global_positive(poly,start,label):
 deg=len(poly)-1;vals=[sum(poly[k]*math.comb(k,j)*start**(k-j)for k in range(j,deg+1))/math.comb(deg,j)for j in range(deg+1)];den=math.lcm(*(v.denominator for v in vals));stack=[([int(v*den)for v in vals],0)];leaves=0
 while stack:
  b,d=stack.pop()
  if min(b)>0:leaves+=1;continue
  assert b[0]>0 and b[-1]>0 and d<100,(label,d)
  l=[b[0]<<deg];r=[b[-1]<<deg]
  for k in range(1,deg+1):b=[x+y for x,y in zip(b,b[1:])];l.append(b[0]<<(deg-k));r.append(b[-1]<<(deg-k))
  for row in (l,r[::-1]):h=math.gcd(*row);stack.append(([x//h for x in row],d+1))
 print(label,'PASS',leaves,flush=True)
global_positive(list(map(F,old['q'])),F(0),'old Fourier factor');global_positive([-x for x in p],tau*E[-1],'new spatial tail')
def ar(x):return arb(x.numerator)/x.denominator
def translated(c,x):return[sum(c[k]*math.comb(k,j)*x**(k-j)for k in range(j,len(c)))for j in range(len(c))]
def radial_bound(left,right):
 order=18;width=right-left;pp=translated(p,tau*left);pp=[x*tau**j for j,x in enumerate(pp)];terms=[(F(1),tau/2,F(0),pp,None),(A,F(240),F(240),None,192),(B,F(120),F(120),None,96)];total=[arb(0)for _ in range(order+1)];rem=arb(0)
 for mul,rate,const,pol,degree in terms:
  if pol is None:pol=[F(math.comb(degree,j))*left**(degree-j)for j in range(degree+1)]
  amp=ar(mul)*ar(const-rate*left).exp();ab=ar(abs(mul))*ar(const-rate*left).exp()
  for k in range(order+1):total[k]+=amp*ar(sum(pol[j]*(-rate)**(k-j)/math.factorial(k-j)for j in range(min(k,len(pol)-1)+1)))
  k=order+1;bound=F(0)
  for j in range(min(k,len(pol)-1)+1):
   sup=F(math.comb(degree,j))*right**(degree-j)if degree is not None else sum(abs(pol[r])*math.comb(r,j)*width**(r-j)for r in range(j,len(pol)))
   bound+=sup*rate**(k-j)/math.factorial(k-j)
  rem+=ab*ar(bound*width**k)
 return[sum((total[j]*ar(width**j*F(math.comb(i,j),math.comb(order,j)))for j in range(i+1)),arb(0))+rem for i in range(order+1)]
leaves=[];cells=0
for i in range(M):
 stack=[('',E[i],E[i+1])];done=[]
 while stack:
  path,l,r=stack.pop();vals=radial_bound(l,r)
  if all(t<=ar(g[i])for t in vals):done.append(path);cells+=1
  else:
   assert len(path)<14,('radial unresolved',i,path)
   mid=(l+r)/2;stack.extend([(path+'0',l,mid),(path+'1',mid,r)])
 leaves.append(done)
 if i%10==0:print('radial',i,'cells',cells,flush=True)
print('all radial bins PASS',M,cells,flush=True)
def capB(t):c=1-2/t;return int(36*(1-c*c)/(1-36*c*c))
def capA(a,t):
 if a+t<=4:return 0
 c=1-2/t;e=(a+t-4)**2/(4*a*t);cc=(c+e)/(1-e)
 return min(capB(t)-1,int(35*(1-cc*cc)/(1-35*cc*cc)))if 35*cc*cc<1 else capB(t)-1
nu=[F(str(x))for x in cand['nu']];assert len(nu)==len(old['angular_cuts']) and min(nu)>=0
weights=[list(map(F,c['weights']))for c in old['angular_cuts']];constants=[F(c['constant'])for c in old['angular_cuts']]
assert all(g[j]>=g[j+1]for j in range(J,M-1));entries=[];budgets=[]
for j in range(-1,J):
 caps=[capB(t)if j<0 else capA(E[j+1],t)for t in E[J+1:]];caps=[min(caps[k:])for k in range(nn)];bd=sum((caps[k]-(caps[k-1]if k else 0))*g[J+k]for k in range(nn))
 if j>=0:bd+=g[j]-sum(nu[k]*weights[k][j]for k in range(len(nu)))
 budgets.append(bd);entries.append(dict(j=j,bound=str(bd)))
kappa=max(budgets);extra=2*sum(v*c for v,c in zip(nu,constants));objective=p[0]/qzero+s+(2*(2**36-1)*kappa+extra)/qzero;density=(tau/8)**18*objective/math.factorial(18)

out=dict(kind='Rational Fourier multiplier plus certified rank-two kernel and class profiles',v=str(v),w=str(w),nu=list(map(str,nu)),extra_objective=str(extra),p=list(map(str,p)),qzero=str(qzero),tau=str(tau),qnew=str(qnew),Z=[[str(s),str(z)],[str(z),str(aZ)]],rho=str(rho),g=list(map(str,g)),E=list(map(str,E)),radial_leaves=leaves,profile_duals=entries,kappa=str(kappa),objective=str(objective),density=str(density),density_decimal=float(density),status='FIRST EXACT REPLAY PASS; independent audit required')
(pth/'certificate.json').write_text(json.dumps(out,indent=2));print('PROFILE and DENSITY PASS',float(density),'objective',float(objective),'cells',cells,flush=True)

assert density<=F('0.000027402464246')
print('INDEPENDENT SAFE ROUNDING PASS',flush=True)
