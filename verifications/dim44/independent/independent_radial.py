"""Independent 44D algebra, centered derivative coupling, and exact rounding.
Uses the frozen rational witness, not discovery output. Global root signs are
also replayed separately by the reviewed complete producer verifier.
"""
import json,math,hashlib,time,sys
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).resolve().parent;S=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else (D if (D/'certificate.json').exists() else D/'frozen_package');raw=(S/'certificate.json').read_bytes();r=json.loads(raw);ang=json.loads((S/'continuum_three_point/cyclic_angular40_exact.json').read_text())
assert r['dimension']==ang['dimension']==44
P=Poly(list(map(Q,r['physical_polynomial'])));H=Poly(list(map(Q,r['fourier_polynomial'])));X=Poly([0,1]);a=Q(22)
# Reconstruct Fourier transform by exact Gaussian-Laguerre eigenbasis elimination.
ls=[Poly([1]),a-2*X]
for k in range(1,P.degree()):ls.append(((2*k+a-2*X)*ls[-1]-(k+a-1)*ls[-2])/Q(k+1))
rem=P;trans=Poly([])
for k in range(P.degree(),-1,-1):
 c=rem[k]/ls[k][k];rem-=c*ls[k];trans+=(-1)**k*c*ls[k]
assert not rem and trans==H and H[0]>0
print('Independent exact eigenbasis Fourier identity PASS',flush=True)
# Angular coefficients by spherical orthogonality rather than triangular extraction.
T=Poly([0,1]);A=Poly(list(map(Q,ang['polynomial_coefficients'])));G=[Poly([1]),T]
for l in range(1,A.degree()):G.append(((2*l+42)*T*G[-1]-l*G[-2])/Q(l+42))
def mean(pol):
 result=Q(0);moment=Q(1)
 for k in range(pol.degree()//2+1):
  if k:moment*=Q(2*k-1,2)/Q(22+k-1)
  result+=pol[2*k]*moment
 return result
cs=[]
for l,g in enumerate(G):
 dim=math.comb(43+l,l)-(math.comb(41+l,l-2)if l>=2 else 0)
 cs.append(dim*mean(A*g))
assert cs==list(map(Q,ang['gegenbauer_coefficients'])) and cs[0]==1 and all(c>=0 for c in cs)
B=sum(cs,Q(0));assert B==A(1)==Q(ang['B_exact'])
print('Independent44D angular orthogonality PSD PASS',flush=True)
s=Q(r['r2']);k=Q(r['kappa']);M=Q(r['M_exact']);assert s==Q(23,4) and k==Q(51,23) and Q(r['R2'])==7 and M>0 and Q(r['b0_exact'])==1
origin=M*B;assert origin*(M/B)==M*M
# Different-length rational Machin enclosure.
def atan_bounds(x,N=96):
 v=sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(N)),Q(0))
 return v,v+x**(2*N+1)/Q(2*N+1)
l1,u1=atan_bounds(Q(1,5));l2,u2=atan_bounds(Q(1,239));pl=16*l1-4*u2;pu=16*u1-4*l2
assert arb(pl)<arb.pi()<arb(pu) and s*k*pl>40
obj=(P[0]+origin)/H[0];assert obj==Q(r['objective_exact'])
upper=Q('1842868599/1000000000000000');assert 0<obj*(pu*s/4)**22/Q(math.factorial(22))<upper
# Independent second-derivative mean-value proof for the normalized-by-h pair expression.
beta=k-1;P1=P.derivative()+beta*P;P2=P.derivative().derivative()+2*beta*P.derivative()+beta*beta*P
AP2=arb_poly([arb(c)for c in P2.coeffs()]);pi=arb.pi();ms=arb(k*s)*pi
stack=[(s*pl,7*pu,0)];leaves=processed=maxdepth=0;start=time.time()
def ev(pol,x):
 v=arb(0)
 for c in reversed(pol.coeffs()):v=v*x+arb(c)
 return v
while stack:
 lo,hi,depth=stack.pop();mid=(lo+hi)/2;half=(hi-lo)/2;m=arb(mid);hh=arb(half)
 ep=(arb(beta)*m-ms).exp();em=(-arb(k)*m+ms).exp()
 f=ev(P,m)*ep-2*arb(M)+arb(M)*em
 fp=ev(P1,m)*ep-arb(M*k)*em
 jets=AP2(arb_poly([m,hh]));p2bound=sum((abs(v)for v in jets),arb(0))
 bound2=p2bound*(arb(beta*hi)-ms).exp()+arb(M*k*k)*(-arb(k*lo)+ms).exp()
 up=f+abs(fp)*hh+bound2*hh*hh/2
 assert up.is_finite() and depth<40 and processed<300000 and time.time()-start<600
 processed+=1;maxdepth=max(maxdepth,depth)
 if up<0:leaves+=1;continue
 assert not f>0,('Positive midpoint',str(mid),str(f))
 stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
 if processed%1000==0:print('coupling',processed,len(stack),depth,time.time()-start,flush=True)
density=arb(obj)*(arb.pi()*arb(s)/4)**22/arb(math.factorial(22));assert density<arb(upper)
out={'source_sha256':hashlib.sha256(raw).hexdigest(),'dimension':44,'independent_eigenbasis_fourier':True,'independent_orthogonality_PSD':True,'origin_rank_one':True,'mean_value_pair_pass':True,'pair_processed':processed,'pair_leaves':leaves,'pair_depth':maxdepth,'coupling_seconds':time.time()-start,'exact_rational_rounding':True,'strict_upper':'0.000001842868599','density_interval':str(density),'scope':'Global radial roots and complete cyclic domains additionally replayed in complete_replay.log'}
(D/'independent_radial.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
