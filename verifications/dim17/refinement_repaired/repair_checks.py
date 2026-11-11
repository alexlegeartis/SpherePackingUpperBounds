"""Independent mixed-Gaussian signs and exact Schur saturation."""
import math
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as Poly,arb

def q(x):return Q(x.numerator,x.denominator)
def chain(p):
 a=[p,p.derivative()]
 while a[-1].degree()>0:
  r=-(a[-2]%a[-1]);assert r
  a.append(r)
 return a

def variation(a,t):
 v=[p[p.degree()] if t is None else p(q(t)) for p in a];s=[int(x>0)-int(x<0)for x in v if x]
 return sum(x!=y for x,y in zip(s,s[1:]))
def no_roots(a,l,h):assert variation(a,l)==variation(a,h)
def bernstein_upper(cs,lo,hi):
 d=len(cs)-1
 powers=[sum(cs[k]*math.comb(k,j)*lo**(k-j)*(hi-lo)**j for k in range(j,d+1))for j in range(d+1)]
 return max(sum(powers[j]*F(math.comb(i,j),math.comb(d,j))for j in range(i+1))for i in range(d+1))
def check_repair(d,ps,hs,As,r2,R2,pl,pu,leaves,depth,pivots):
 g=d['slow_gaussian_correction'];eta=F(g['eta']);beta=F(g['beta']);c=F(g['c']);n=d['dimension'];a=F(n,2)
 assert n==17 and beta==F(1,4) and eta==F(13,10**15) and c==F(667,100)
 assert len(ps)==len(hs)==32 and pl*r2>c
 # Fourier of (c-X)exp(-beta X) is beta^(-n/2)
 # times (c-n/(2 beta)+Y/beta²)exp(-Y/beta).
 factor=F(2**n);cutoff=a*beta-c*beta**2;cost=eta*factor*(a/beta-c)
 assert cutoff==F(2733,1600) and cost>0
 H=Poly(list(map(q,hs)));P=Poly(list(map(q,ps)))
 h=H-q(cost);hc=chain(h);assert h[0]>0 and h(q(cutoff))>0;no_roots(hc,F(0),cutoff)
 # On [0,cutoff], corrected Fourier /exp(-Y)>=H-cost;
 # beyond cutoff, additive Fourier term is nonnegative.
 pc=chain(P);left=R2*pl;assert left>c;reports=[]
 for lo,hi in [tuple(map(F,z))for z in g['repair_intervals']]:
  assert left<lo<hi and P(q(left))<0 and P(q(lo))<0 and P(q(hi))<0
  no_roots(pc,left,lo)
  upper=bernstein_upper(ps,lo,hi)
  allowance=arb(q(eta*(lo-c)))*(arb(q((1-beta)*lo))).exp()
  assert arb(q(upper))<allowance
  reports.append({'interval':[str(lo),str(hi)],'polynomial_upper':str(arb(q(upper))),'allowance_lower':str(allowance)})
  left=hi
 no_roots(pc,left,None);assert P(q(left))<0
 # The physical correction is strictly negative at all allowed separations.
 norm=hs[0]-cost;assert norm>0 and norm==F(g['corrected_normalization'])
 numerator=ps[0]+eta*c+As[0][0][0];obj=numerator/norm;assert obj==F(d['exact_objective'])
 vol=arb.pi()**8*arb(q(r2))**8*arb(q(r2)).sqrt()*2*math.factorial(9)/math.factorial(18)
 density=arb(q(obj))*vol
 sq=lambda o:o**2*pu**16*r2**17*4*math.factorial(9)**2/math.factorial(18)**2
 assert obj>0 and sq(obj)<F('.018250915')**2
 # Solve B x=b exactly and reconstruct the saturated PSD matrix W B W^T.
 A=As[0];B=[row[1:]for row in A[1:]];b=[row[0]for row in A[1:]];m=len(b);aug=[B[i][:]+[b[i]]for i in range(m)]
 for j in range(m):
  pivot=aug[j][j];assert pivot>0
  for i in range(j+1,m):
   z=aug[i][j]/pivot
   for k in range(j,m+1):aug[i][k]-=z*aug[j][k]
 x=[F(0)]*m
 for i in range(m-1,-1,-1):x[i]=(aug[i][-1]-sum(aug[i][j]*x[j]for j in range(i+1,m)))/aug[i][i]
 assert all(sum(B[i][j]*x[j]for j in range(m))==b[i]for i in range(m))
 alpha=A[0][0]-sum(b[i]*x[i]for i in range(m));assert alpha>0
 W=[x]+[[F(i==j)for j in range(m)]for i in range(m)]
 assert all(A[i][j]-(alpha if i==j==0 else 0)==sum(W[i][k]*B[k][l]*W[j][l]for k in range(m)for l in range(m))for i in range(m+1)for j in range(m+1))
 sharp=(numerator-alpha)/norm;assert sharp>0 and sq(sharp)<F('.018204096')**2
 return {'pass':True,'positive_pivots':pivots,'pair_leaves':leaves,'pair_depth':depth,'repair_intervals':reports,'density_interval':str(density),'strict_upper':'.018250915','corrected_normalization':str(norm),'alpha_exact':str(alpha),'schur_density_interval':str(arb(q(sharp))*vol),'schur_strict_upper':'.018204096'}
