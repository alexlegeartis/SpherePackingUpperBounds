"""Exact initial Bernstein coefficients; outward-rounded dyadic de Casteljau.
Checks K(r,s,t)<=0 for 1<=r,s<=R and distance>=R.
"""
import sys,json,time,math,numpy as np
sys.path.insert(0,'solver_packages')
from fractions import Fraction
from flint import fmpq,fmpq_poly,fmpq_mat,arb,ctx
ctx.prec=128
path=sys.argv[1];data=json.load(open(path));src=json.load(open(data['source']));p=src['p'];K=src['K'];R=fmpq(str(Fraction(str(src['R']))));Qs=[fmpq_mat([[fmpq(z)for z in row]for row in M])for M in data['matrices']]
assert src["J"]==1
X=fmpq_poly([0,1]);Z=2*X-1;Ts=[fmpq_poly([1]),Z]
for i in range(1,p):Ts.append(2*Z*Ts[-1]-Ts[-2])
top=fmpq(1,2);T=-1+(top+1)*X;Us=[fmpq_poly([1]),2*T]
for k in range(1,K):Us.append(2*T*Us[-1]-Us[-2])
def bern(poly,d):return [sum((poly[j]*fmpq(math.comb(i,j),math.comb(d,j))for j in range(min(i,poly.degree())+1)),fmpq(0))for i in range(d+1)]
Br=fmpq_mat([bern(Ts[i],p)for i in range(p+1)]);Ba=[bern(Us[k]/(k+1),K)for k in range(K+1)]
M=[]
for k,Q in enumerate(Qs):
 q=fmpq_mat([[Q[i+(k==0),j+(k==0)]for j in range(p+1)]for i in range(p+1)]);M.append(Br.transpose()*q*Br)
lo=np.empty((p+1,p+1,K+1));hi=lo.copy()
for i in range(p+1):
 for j in range(p+1):
  for l in range(K+1):
   z=sum((M[k][i,j]*Ba[k][l]for k in range(K+1)),fmpq(0));v=arb(z)
   lo[i,j,l]=np.nextafter(float(v.lower()),-np.inf);hi[i,j,l]=np.nextafter(float(v.upper()),np.inf)
print('INITIAL',float(np.min(lo)),float(np.max(hi)),flush=True)
def split(l,h,axis):
 a=np.moveaxis(l,axis,0).copy();b=np.moveaxis(h,axis,0).copy();d=a.shape[0]-1
 ll=np.empty_like(a);lh=np.empty_like(b);rl=np.empty_like(a);rh=np.empty_like(b)
 ll[0]=a[0];lh[0]=b[0];rl[d]=a[-1];rh[d]=b[-1]
 for k in range(1,d+1):
  a=np.nextafter(np.nextafter(a[:-1]+a[1:],-np.inf)*.5,-np.inf)
  b=np.nextafter(np.nextafter(b[:-1]+b[1:],np.inf)*.5,np.inf)
  ll[k]=a[0];lh[k]=b[0];rl[d-k]=a[-1];rh[d-k]=b[-1]
 return tuple(np.moveaxis(z,0,axis)for z in [ll,lh,rl,rh])
Ra=arb(R);R2=arb(R*R);ta=arb(top+1)
def outside(box):
 # Distance is decreasing in t and convex separately in r and s.
 # Its maximum on each radial rectangle is at a corner, with lowest t.
 t=-1+ta*arb(box[2][0])
 for u in box[0]:
  r=1+(Ra-1)*arb(u)
  for v in box[1]:
   s=1+(Ra-1)*arb(v)
   if not (r*r+s*s-2*r*s*t<R2):return False
 return True
stack=[(lo,hi,((0.,1.),)*3,0)];checked=passed=excluded=0;maxdepth=0;fail=None;t0=time.time();limit=int(sys.argv[2])if len(sys.argv)>2 else 500000
while stack:
 l,h,box,d=stack.pop();checked+=1;maxdepth=max(maxdepth,d)
 if np.max(h)<0:passed+=1;continue
 if outside(box):excluded+=1;continue
 if checked>=limit or d>=55:
  fail={'box':box,'depth':d,'upper':float(np.max(h)),'remaining':len(stack)};break
 mid=(l+h)*.5;scores=[float(np.max(np.abs(np.diff(mid,axis=i))))for i in range(3)];axis=int(np.argmax(scores))
 ll,lh,rl,rh=split(l,h,axis);a,b=box[axis];m=(a+b)*.5;bx1=list(box);bx2=list(box);bx1[axis]=(a,m);bx2[axis]=(m,b)
 stack.append((rl,rh,tuple(bx2),d+1));stack.append((ll,lh,tuple(bx1),d+1))
 if checked%10000==0:print('PROGRESS',checked,'stack',len(stack),'seconds',time.time()-t0,flush=True)
out={'angular_top':'1/2','verified':fail is None and not stack,'checked':checked,'passed':passed,'excluded':excluded,'max_depth':maxdepth,'seconds':time.time()-t0,'failure':fail,'source':path,'method':'Exact rational initial tensor Bernstein coefficients, outward-rounded float64 dyadic subdivision, Arb domain exclusion'}
print(json.dumps(out),flush=True);json.dump(out,open(path.replace('_psd_exact.json','_far_bernstein.json'),'w'),indent=2)
