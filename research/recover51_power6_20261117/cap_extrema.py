"""Root-isolated extrema for discovery cuts; not a bound certificate."""
from fractions import Fraction
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
import itertools

def extrema(co,B,n=51):
 ctx.prec=384;co=list(map(Q,co));B=[[[Q(str(Fraction(str(x))))for x in row]for row in b]for b in B];x=P([0,1]);G=[P([1]),x]
 for k in range(2,len(co)):G.append(((2*k+n-4)*x*G[-1]-(k-1)*G[-2])/(k+n-3))
 T=[[Q(1),Q(0),Q(0)],[Q(1),Q(1,2),Q(0)],[Q(1),Q(1),Q(1)]];points=[];maximum=arb(-10**100)
 for a in range(3):
  for b in range(a,3):
   poly=sum((co[k]*sum((T[a][i]*B[k][i][j]*T[b][j]for i in range(3)for j in range(3)),Q(0))*G[k]for k in range(len(co))),P([0]));ap=arb_poly([arb(v)for v in poly]);roots=poly.derivative().complex_roots();assert sum(m for z,m in roots)==poly.degree()-1
   maximum=max(maximum,arb(poly(-1)).upper(),arb(poly(Q(1,2))).upper())
   for z,m in roots:
    if z.imag.contains(0) and not z.real < -1 and not z.real>arb(1)/2:points.append(float(z.real.mid()));maximum=max(maximum,ap(z.real).upper())
 return points,float(maximum)
