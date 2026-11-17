import sys,json,math
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_poly as P
if not __debug__:raise RuntimeError('Assertions required')
D=Path(__file__).resolve().parent;r=json.loads((D/'cap_eps1000.json').read_text());n=r['dimension'];co=[Q(str(Fraction(str(v))))for v in r['coefficients']];assert co[0]==1 and min(co)>=0;s=Q(r['s']);assert -1<s<1
x=P([0,1]);g=[P([1]),x]
for k in range(2,len(co)):g.append(((2*k+n-4)*x*g[-1]-(k-1)*g[-2])/(k+n-3))
f=sum((c*gp for c,gp in zip(co,g)),P([0]));deg=f.degree();stack=[(Q(-1),s,0)];nodes=0;leaves=0
while stack:
 lo,hi,depth=stack.pop();nodes+=1;p=f(P([lo,hi-lo]));bs=[sum((p[j]*Q(math.comb(i,j),math.comb(deg,j))for j in range(i+1)),Q(0))for i in range(deg+1)]
 if max(bs)<=0:leaves+=1;continue
 assert depth<45,(lo,hi,max(bs))
 mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
assert nodes==2*leaves-1;B=sum(co);assert f(1)==B;C=int(B.numerator//B.denominator)
out=dict(status='EXACT SPHERICAL CODE CAP PASS',dimension=n,s=str(s),epsilon=r['epsilon'],cap=C,ratio=str(B),gegenbauer=list(map(str,co)),polynomial=list(map(str,f)),nodes=nodes,leaves=leaves);(D/'contact_cap_certificate.json').write_text(json.dumps(out,indent=2));print(n,'CAP PASS',C,nodes,flush=True)
