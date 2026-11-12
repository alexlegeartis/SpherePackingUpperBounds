"""Deterministic exact export of frozen reconstruction inputs; no optimization."""
import json,sys,hashlib,argparse
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=D/'exported');args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
def rat(x):return Q(str(F(float(x))))
def save(name,j):(out/name).write_text(json.dumps(j,indent=2))
a=json.loads((D/'inputs/angular_discovery.json').read_text());c=list(map(rat,a['coefficients']));c[0]-=Q(1,10000);c=[z/c[0] for z in c];x=P([0,1]);gs=[P([1]),x]
for k in range(1,len(c)-1):gs.append(((2*k+8)*x*gs[-1]-k*gs[-2])/(k+8))
poly=sum((v*g for v,g in zip(c,gs)),P());angular=dict(dimension=10,a_exact='8',b_exact='16',gegenbauer_coefficients=list(map(str,c)),B_exact=str(sum(c)),polynomial_coefficients=list(map(str,poly.coeffs())))
save('angular.json',angular)
r=json.loads((D/'inputs/refit10.json').read_text());gs=[P([1]),5-2*x]
for k in range(1,len(r['radial'])-1):gs.append(((2*k+5-2*x)*gs[-1]-(k+4)*gs[-2])/(k+1))
gs=[g/g(0) for g in gs];p=sum((rat(v)*g for v,g in zip(r['radial'],gs)),P());h=sum(((-1)**k*rat(v)*g for k,(v,g) in enumerate(zip(r['radial'],gs))),P());N=h(0);p/=N;h/=N;C=rat(r['C'])/N;eps=Q(1,10**6)
p=(p+eps*(1-2*x/5))/(1-eps);h=(h-eps*(1-2*x/5))/(1-eps);C=C/(1-eps)*Q(1000001,1000000)
cert=dict(dimension=10,physical=list(map(str,p.coeffs())),fourier=list(map(str,h.coeffs())),C_exact=str(C),epsilon_exact=str(eps),s_exact='5/4')
save('certificate.json',cert);(out/'fit.json').write_bytes((D/'inputs/fit10.json').read_bytes())
print('Exported exact proof inputs to',out)
