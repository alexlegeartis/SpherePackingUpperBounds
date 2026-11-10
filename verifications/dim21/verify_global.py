import sys,json,math,time
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_CEILING
from flint import fmpq,fmpq_poly,arb,ctx
ctx.prec=256
r=json.load(open(sys.argv[1] if len(sys.argv)>1 else Path(__file__).with_name('certificate.json')))
n=int(r['dimension']);K=int(r['radial_nodes']);L=int(r['angular_degree']);D=int(r['radial_degree']);r2=F(r['exclusion_squared']);R2=F(r['truncation_squared'])
assert n==21 and L==16 and K==5 and r2==F(21,8) and R2==F(63,20) and R2<2*r2
Q=lambda v:fmpq(str(v))
As=[[[F(z)for z in row]for row in A]for A in r['matrices']]
X=fmpq_poly([0,1]);P=fmpq_poly(list(map(fmpq,r['physical_polynomial'])));H=fmpq_poly(list(map(fmpq,r['fourier_polynomial'])))
def ldlt(A):
 k=len(A);d=[];ll=[[F(int(i==j))for j in range(k)]for i in range(k)]
 for j in range(k):
  dj=A[j][j]-sum(ll[j][s]**2*d[s]for s in range(j));d.append(dj)
  if dj<=0:return False,float(dj)
  for i in range(j+1,k):ll[i][j]=(A[i][j]-sum(ll[i][s]*ll[j][s]*d[s]for s in range(j)))/dj
 return True,min(map(float,d))

assert all(ldlt(A)[0]for A in As)
print('All 17 rational PSD blocks: PASS',flush=True)
m=fmpq(n,2);ls=[fmpq_poly([1]),m-2*X]
for k in range(1,D):ls.append(((2*k+m-2*X)*ls[-1]-(k+m-1)*ls[-2])/fmpq(k+1))
rem=P;fh=fmpq_poly([])
for k in range(D,-1,-1):
 c=rem[k]/ls[k][k];rem-=c*ls[k];fh+=((-1)**k)*c*ls[k]
assert not rem and fh==H and H[0]==fmpq(r['normalization'])
print('Exact Fourier transform identity: PASS',flush=True)
def sturm(p):
 seq=[p,p.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);num=rem.numer();cont=abs(num.content());seq.append(fmpq_poly(num/cont))
 return seq
def variation(vals):
 ss=[1 if a>0 else -1 for a in vals if a];return sum(a!=b for a,b in zip(ss,ss[1:]))
def roots_pos(p):
 seq=sturm(p);return variation([p[0]for p in seq])-variation([p[p.degree()]for p in seq])
pi_lo=F(314159265358979323846264338327950288419716939937510,10**50)
assert roots_pos(H)==0 and H[0]>0
print('Fourier sign: PASS',flush=True)
assert arb(Q(pi_lo))<arb.pi()
assert roots_pos(-P(X+Q(R2*pi_lo)))==0 and P(Q(R2*pi_lo))<0
print('Exterior physical sign: PASS',flush=True)
# Angular rectangle certificate
G=int(r['subdivisions_per_radial_cell']);nodes=[r2+(R2-r2)*i/(G*(K-1))for i in range(G*(K-1)+1)]
ps=[fmpq_poly([1]),X]
for l in range(1,L):ps.append(((2*l+n-2)*X*ps[-1]-l*ps[-2])/fmpq(l+n-2))
def w(u):
 q=(u-r2)/(R2-r2)*(K-1);i=min(K-2,max(0,int(q)));v=[F(0)]*K;v[i]=1-(q-i);v[i+1]=q-i;return v
def bil(a,A,b):return sum(a[i]*A[i][j]*b[j]for i in range(K)for j in range(K))
def cup(u,v):
 with localcontext()as c:
  c.prec=80;dd=lambda z:Decimal(z.numerator)/Decimal(z.denominator);z=dd(u+v-r2)/(2*dd(u*v).sqrt());ans=F(int((z*Decimal(10**60)).to_integral_value(rounding=ROUND_CEILING)),10**60)
  while 4*ans*ans*u*v<(u+v-r2)**2:ans+=F(1,10**60)
  assert ans>0 and 4*ans*ans*u*v>=(u+v-r2)**2
  return ans
failed=[];checked=0;t0=time.time()
for i,u in enumerate(nodes):
 for j in range(i,len(nodes)):
  v=nodes[j];a=w(u);b=w(v);cs=[bil(a,[rr[1:]for rr in As[0][1:]],b)]+[bil(a,A,b)for A in As[1:]];poly=sum((Q(c)*p for c,p in zip(cs,ps)),start=fmpq_poly([]));cu=cup(nodes[min(i+1,len(nodes)-1)],nodes[min(j+1,len(nodes)-1)]);seq=sturm(poly);num=variation([p(-1)for p in seq])-variation([p(Q(cu))for p in seq]);ok=num==0 and poly(-1)<0 and poly(Q(cu))<0
  checked+=1
  if not ok:failed.append([i,j,num,str(poly(-1)),str(poly(Q(cu)))])
print('EXACT_ANGULAR',checked,'failures',len(failed),'seconds',time.time()-t0,flush=True)

assert not failed and checked==561
objective=(P[0]+Q(As[0][0][0]))/H[0]
assert objective==fmpq(r['exact_objective'])
density=arb.pi()**10*math.factorial(11)*4**11/math.factorial(22)*(arb(fmpq(21,32))**21).sqrt()*arb(objective)
assert density<arb(r['certified_density_upper'])
print('Density enclosure:',density,flush=True)
print('Global algebraic and angular checks: PASS',flush=True)
