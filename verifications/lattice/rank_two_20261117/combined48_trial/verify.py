import json,math,sys,hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
ctx.prec=768;sys.set_int_max_str_digits(0)
pth=Path(__file__).parent;old=json.loads((pth.parent/'certificate48.json').read_text());cand=json.loads((pth/'input_candidate.json').read_text());basis=json.loads((pth/'input_basis.json').read_text());tau=F(basis['tau']);qzero=F(basis['qzero']);qnew=F(basis['q']);polys=[list(map(F,v))for v in basis['polys']];assert tau==F(old['tau']) and qzero==F(old['q'][0]);assert polys[0]==list(map(F,old['p']))
co=[F(str(v))for v in json.loads((pth.parent/'rank2_candidate_B783489561_20261117.json').read_text())['even_coefficients']];co[0]-=F(2,10000);assert qnew==sum(co)/co[0]
def op(p):return [((48+4*j)*(p[j]if j<len(p)else 0)-(p[j-1]if j>=1 else 0)-(j+1)*(96+4*j)*(p[j+1]if j+1<len(p)else 0))/tau for j in range(len(p)+1)]
assert op(polys[0])==polys[1] and op(polys[1])==polys[2]
import sympy as sp
def matrix(key):
 raw=[[F(str(v))for v in row]for row in cand[key]];n=len(raw)
 return [[(raw[i][j]+raw[j][i])/2 for j in range(n)]for i in range(n)]
def psd(M):
 assert all(sp.Matrix(M)[:k,:k].det()>0 for k in range(1,len(M)+1))
FM=matrix('fourier_psd');FN=matrix('fourier_halfline_psd');VM=matrix('fourier_multiplier_psd');psd(FM);psd(FN)
m=[F(0)]*9
for mat,shift in [(FM,0),(FN,1)]:
 for i,row in enumerate(mat):
  for j,c in enumerate(row):
   for k in range(i+j+1):m[k+shift]+=c*math.comb(i+j,k)*(-1)**(i+j-k)
VM[0][0]=1-m[0];psd(VM);v=2*VM[0][1];w=VM[1][1];coeffs=[VM[0][0],v,w]+m
assert VM[0][0]+m[0]==1
expected=[F(1)]
for j in range(9):
 assert polys[j+3]==[qzero*x for x in expected]
 expected=op(expected)
p=[sum(z*(poly[j]if j<len(poly)else 0)for z,poly in zip(coeffs,polys))for j in range(max(map(len,polys)))];K=len(p)-1

Z=[[F(str(z))for z in row]for row in cand['origin']];assert Z[0][1]==Z[1][0];z=Z[0][1]-F(1,10**8);aZ=Z[1][1];assert aZ>0 and z<0;s=z*z/aZ+F(1,10**20);assert s*aZ-z*z>0
rho=F(1,10**70);assert F(168**168,math.factorial(168))>1/rho;factor=(1+rho)/(1-rho);A=2*qzero*aZ*factor/qnew;B=2*qzero*z/qnew;assert A>0 and B<0 and A+B<0
E=list(map(F,old['E']));g=[max(F(0),F(str(x)))+F(1,10**24)for x in cand['g']];M=len(g);J=E.index(F(2));nn=M-J;assert len(E)==M+1 and not old['angular_cuts'];assert len(cand['profile_y'])==J+1
ds=list(map(F,old['d']));oldp=list(map(F,old['p']));oldq=list(map(F,old['q']));N=len(ds)-1
for j in range(N+1):
 den=math.factorial(j)*math.prod(F(24)+a for a in range(j));values=[ds[k]*F(math.factorial(k),math.factorial(k-j))for k in range(j,N+1)]
 assert oldp[j]==(-1)**j*sum(values)/den
 assert oldq[j]==(-1)**j*sum(a*(-1)**k for k,a in zip(range(j,N+1),values))/den
assert hashlib.sha256((pth.parent/'rank2_candidate_B783489561_20261117.json').read_bytes()).hexdigest()=='119dc3c5ef12da19e6edb01332e87e4c472da94d6c60c23fa438f74b7da147ae'
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
def capB(t):c=1-2/t;return int(48*(1-c*c)/(1-48*c*c))
def capA(a,t):
 if a+t<=4:return 0
 c=1-2/t;e=(a+t-4)**2/(4*a*t);cc=(c+e)/(1-e)
 return min(capB(t)-1,int(47*(1-cc*cc)/(1-47*cc*cc)))if 47*cc*cc<1 else capB(t)-1
entries=[];budgets=[]
for idx,yy in enumerate(cand['profile_y']):
 j=idx-1;caps=[capB(t)if j<0 else capA(E[j+1],t)for t in E[J+1:]];caps=[min(caps[k:])for k in range(nn)];rhs=list(map(F,caps));cols=[[F(int(k<=h))for h in range(nn)]for k in range(nn)]
 for fac in ([F(1),F(21,20)]if j>=0 else[F(1)]):
  for h,t in enumerate(E[J+1:]):
   alpha=fac*47*(t-2)/(4*(t-1));condition=j>=0 and(fac>1 or E[j]>=F(191,96));rhs.append(47+alpha*(E[j+1]-2)if condition else F(48))
   for k in range(nn):cols[k].append(1-alpha*(E[J+k+1]-2)if k<=h else F(0))
 y=[max(F(0),F(str(z)))for z in yy];assert len(y)==len(rhs);deficit=max([F(0)]+[g[J+k]-sum(a*b for a,b in zip(y,cols[k]))for k in range(nn)]);y[nn-1]+=deficit+F(1,10**30)
 assert all(sum(a*b for a,b in zip(y,cols[k]))>=g[J+k]for k in range(nn));bd=sum(a*b for a,b in zip(y,rhs))+(g[j]if j>=0 else F(0));budgets.append(bd);entries.append(dict(j=j,y=list(map(str,y)),bound=str(bd)))
kappa=max(budgets);objective=p[0]/qzero+s+2*(2**48-1)*kappa/qzero;density=(tau/8)**24*objective/math.factorial(24)
out=dict(kind='Rational Fourier multiplier plus certified rank-two kernel and class profiles',v=str(v),w=str(w),fourier_multiplier=list(map(str,coeffs[:3])),fourier_sos=[[str(v)for v in row]for row in FM],fourier_halfline_sos=[[str(v)for v in row]for row in FN],fourier_additive_coefficients=list(map(str,m)),p=list(map(str,p)),qzero=str(qzero),tau=str(tau),qnew=str(qnew),Z=[[str(s),str(z)],[str(z),str(aZ)]],rho=str(rho),g=list(map(str,g)),E=list(map(str,E)),radial_leaves=leaves,profile_duals=entries,kappa=str(kappa),objective=str(objective),density=str(density),density_decimal=float(density),status='FIRST EXACT REPLAY PASS; independent audit required')
(pth/'certificate.json').write_text(json.dumps(out,indent=2));print('PROFILE and DENSITY PASS',float(density),'objective',float(objective),'cells',cells,flush=True)
