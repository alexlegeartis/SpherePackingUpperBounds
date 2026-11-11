import sys,json,time,math
sys.path.insert(0,'solver_packages')
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
ctx.prec=128
j=json.load(open('continuum_three_point/cyclic_angular40_exact.json'));Q=fmpq_poly(list(map(fmpq,j['polynomial_coefficients'])));D=Q.degree();B=arb(fmpq(j['B_exact']));AQ=arb_poly([arb(x)for x in Q.coeffs()]);alpha=int(sys.argv[1])if len(sys.argv)>1 else 64
def bern(a,b):
 p=Q(fmpq_poly([a,b-a]));return tuple(sum((p[k]*fmpq(math.comb(i,k),math.comb(D,k))for k in range(i+1)),fmpq(0))for i in range(D+1))
def split(v):
 left=[v[0]];right=[v[-1]]
 while len(v)>1:
  v=tuple((v[i]+v[i+1])/2 for i in range(len(v)-1));left.append(v[0]);right.append(v[-1])
 return tuple(left),tuple(reversed(right))
cache={}
def get(axis,i,level):
 key=axis,i,level
 if key in cache:return cache[key]
 if level==0:
  a,b=(fmpq(0),fmpq(1,2))if axis==0 else(fmpq(501,1000),fmpq(1));r=(a,b,bern(a,b))
 else:
  a,b,v=get(axis,i//2,level-1);m=(a+b)/2;l,r=split(v);cache[axis,2*(i//2),level]=(a,m,l);cache[axis,2*(i//2)+1,level]=(m,b,r);return cache[key]
 cache[key]=r;return r

def qr(x):
 a=max(x.lower(),arb(-1));b=min(x.upper(),arb(1));c=(a+b)/2;r=(b-a)/2;T=AQ(arb_poly([c,r]));err=sum((abs(T[i])for i in range(1,len(T))),arb(0));return (T[0]-err).union(T[0]+err)
stack=[(0,0,0,0)];passed=0;maxlev=[0,0];start=time.time();skipped=0
while stack:
 i,l,k,m=stack.pop();ta,tb,T=get(0,i,l);ua=fmpq(1)+fmpq(4*k,2**m);ub=fmpq(1)+fmpq(4*(k+1),2**m)
 if 2*ta*ua>1:skipped+=1;continue
 U=max(T);t=arb(ta).union(arb(tb));u=arb(ua).union(arb(ub));v=(1+u*u-2*u*t).sqrt()
 if not v.is_finite():val=arb('nan')
 else:
  z1=(1-u*t)/v;z2=(u-t)/v;g1=(1-2*u*t);g2=(u*u-2*u*t);val=arb(U)+qr(z1)*(-alpha*g1).exp()+qr(z2)*(-alpha*g2).exp()
 if val<0:passed+=1;maxlev=[max(l,maxlev[0]),max(m,maxlev[1])];continue
 if max(l,m)>=40:raise RuntimeError((i,l,k,m,str(val)))
 if U>=0 or l<=m:stack.extend([(2*i,l+1,k,m),(2*i+1,l+1,k,m)])
 else:stack.extend([(i,l,2*k,m+1),(i,l,2*k+1,m+1)])
 if passed and passed%10000==0:print(passed,len(stack),maxlev,time.time()-start,flush=True)
print('PASS',alpha,passed,maxlev,time.time()-start,flush=True)
json.dump({'pass':True,'alpha':alpha,'rectangles':passed,'skipped':skipped,'max_levels':maxlev,'precision':ctx.prec,'elapsed':time.time()-start,'scope':'Finite t,u triangle domain only; use separate analytic t<=0/u>=5 tail.'},open('continuum_three_point/cyclic_optimized_actual'+str(alpha)+'_exact.json','w'),indent=2)
