"""Exact polynomial numerator and tensor-Bernstein continuum certificate."""
if not __debug__: raise RuntimeError("Verification must not run under python -O.")

from pathlib import Path
import sys,json,time,math,numpy as np
from fractions import Fraction

from flint import fmpq,fmpq_poly,fmpq_mat,arb,ctx
ctx.prec=128;path=sys.argv[1];data=json.load(open(path));src=json.load(open((Path(path).resolve().parent/data['source']) if not Path(data['source']).is_absolute() else data['source']));p=src['p'];K=src['K'];R=fmpq(str(Fraction(str(src['R']))));assert src['J']==1 and R<2
n=src["dimension"]
Qs=[fmpq_mat([[fmpq(z)for z in row]for row in M])for M in data['matrices']]
X=fmpq_poly([0,1]);Z=(2*X-R-1)/(R-1);Ts=[fmpq_poly([1]),Z]
for i in range(1,p):Ts.append(2*Z*Ts[-1]-Ts[-2])
G=[fmpq_poly([1]),X]
for k in range(1,K):G.append(((2*k+n-2)*X*G[-1]-k*G[-2])/fmpq(k+n-2))
B=fmpq_mat([[Ts[i][j]for j in range(p+1)]for i in range(p+1)]);D={};t0=time.time()
for k,Q in enumerate(Qs):
 off=(k==0);qr=fmpq_mat([[Q[i+off,j+off]for j in range(p+1)]for i in range(p+1)]);M=B.transpose()*qr*B
 for j in range(k//2+1):
  m=k-2*j;cc=G[k][m]/(2**m)
  for u in range(m+1):
   for v in range(m-u+1):
    w=m-u-v;c=cc*((-1)**w)*math.comb(m,u)*math.comb(m-u,v)
    for a in range(p+1):
     for b in range(p+1):
      z=c*M[a,b]
      if not z:continue
      e=(K-m+a+2*u,K-m+b+2*v,K+2*w)
      for q in [e,(e[0],e[2],e[1]),(e[2],e[0],e[1])]:D[q]=D.get(q,fmpq(0))+z
print('POLYNOMIAL',len(D),'seconds',time.time()-t0,flush=True)
# Independent exact evaluation at a nonsymmetric rational triangle.
def kval(a,b,c):
 v=[t(a)for t in Ts];w=[t(b)for t in Ts];t=(a*a+b*b-c*c)/(2*a*b);U=[fmpq(1),t]
 for k in range(1,K):U.append(((2*k+n-2)*t*U[-1]-k*U[-2])/fmpq(k+n-2))
 ans=fmpq(0)
 for k,Q in enumerate(Qs):
  o=(k==0);q=sum((v[i]*Q[i+o,j+o]*w[j]for i in range(p+1)for j in range(p+1)),fmpq(0));ans+=q*U[k]
 return ans
a,b,c=fmpq(1),R,(1+R)/2;pw=[[z**i for i in range(3*K+p+1)]for z in [a,b,c]];ev=sum((z*pw[0][i]*pw[1][j]*pw[2][k]for(i,j,k),z in D.items()),fmpq(0));direct=(kval(a,b,c)+kval(a,c,b)+kval(b,c,a))*(a*b*c)**K
assert ev==direct,'polynomial construction mismatch';print('EXACT_IDENTITY_TEST_PASS',flush=True)
d=max(max(e)for e in D);n=d+1;delta=R-1
BT=fmpq_mat([[sum((fmpq(math.comb(e,j)*math.comb(i,j),math.comb(d,j))*delta**j for j in range(min(i,e)+1)),fmpq(0))for e in range(n)]for i in range(n)])
A=np.empty((n,n,n),dtype=object);A.fill(fmpq(0))
for e,v in D.items():A[e]=v
for axis in range(3):
 V=np.moveaxis(A,axis,0).reshape(n,n*n);M=BT*fmpq_mat(V.tolist());V=np.array([[M[i,j]for j in range(n*n)]for i in range(n)],dtype=object);A=np.moveaxis(V.reshape(n,n,n),0,axis)
 print('TRANSFORM',axis,'seconds',time.time()-t0,flush=True)
lo=np.empty(A.shape);hi=np.empty(A.shape)
for e in np.ndindex(A.shape):
 z=arb(A[e]);lo[e]=np.nextafter(float(z.lower()),-np.inf);hi[e]=np.nextafter(float(z.upper()),np.inf)
print('INITIAL',A.shape,float(np.min(lo)),float(np.max(hi)),flush=True);del A,D
# Save initialized rigorous interval enclosure to permit independent reruns.
np.savez_compressed(path.replace('_psd_exact.json','_cyclic_bernstein_initial.npz'),lo=lo,hi=hi)
def split(l,h,axis):
 a=np.moveaxis(l,axis,0).copy();b=np.moveaxis(h,axis,0).copy();d=a.shape[0]-1;ll=np.empty_like(a);lh=np.empty_like(b);rl=np.empty_like(a);rh=np.empty_like(b);ll[0]=a[0];lh[0]=b[0];rl[d]=a[-1];rh[d]=b[-1]
 for k in range(1,d+1):
  a=np.nextafter(np.nextafter(a[:-1]+a[1:],-np.inf)*.5,-np.inf);b=np.nextafter(np.nextafter(b[:-1]+b[1:],np.inf)*.5,np.inf)
  ll[k]=a[0];lh[k]=b[0];rl[d-k]=a[-1];rh[d-k]=b[-1]
 return tuple(np.moveaxis(z,0,axis)for z in [ll,lh,rl,rh])
stack=[(lo,hi,((Fraction(0),Fraction(1)),)*3,0)];checked=passed=excluded=0;maxdepth=0;fail=None;limit=int(sys.argv[2])if len(sys.argv)>2 else 50000
while stack:
 l,h,box,depth=stack.pop();checked+=1;maxdepth=max(maxdepth,depth)
 if np.max(h)<0:passed+=1;continue
 if box[0][0]>box[1][1]or box[1][0]>box[2][1]:excluded+=1;continue
 if checked>=limit or depth>=65:
  fail={'box':box,'depth':depth,'upper':float(np.max(h)),'remaining':len(stack)};break
 mid=(l+h)*.5;scores=[float(np.max(np.abs(np.diff(mid,axis=i))))for i in range(3)];axis=int(np.argmax(scores));ll,lh,rl,rh=split(l,h,axis);a,b=box[axis];m=(a+b)/2;bx1=list(box);bx2=list(box);bx1[axis]=(a,m);bx2[axis]=(m,b);stack.extend([(rl,rh,tuple(bx2),depth+1),(ll,lh,tuple(bx1),depth+1)])
 if checked%100==0:print('PROGRESS',checked,'stack',len(stack),'seconds',time.time()-t0,flush=True)
out={'verified':fail is None and not stack,'checked':checked,'passed':passed,'excluded':excluded,'max_depth':maxdepth,'seconds':time.time()-t0,'failure':fail,'source':path,'method':'Exact rational numerator (abc)^K times cyclic sum; exact rational Bernstein initialization; outward-rounded dyadic subdivision; permutation symmetry'}
print(json.dumps(out,default=str),flush=True);json.dump(out,open(path.replace('_psd_exact.json','_cyclic_bernstein.json'),'w'),indent=2,default=str)
