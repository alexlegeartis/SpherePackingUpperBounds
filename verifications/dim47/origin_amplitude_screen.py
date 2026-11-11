from pathlib import Path
from fractions import Fraction as F
import mpmath as mp
import numpy as np,json
from scipy.optimize import minimize_scalar
D=Path(__file__).parent;S=D/'continuum_three_point';mp.mp.dps=100
r=json.loads((S/'adaptive_cyclic40_radial111.json').read_text());a=json.loads((S/'cyclic_angular40_exact.json').read_text());B=F(a['B_exact']);M=F(r['neighbor_cost'])/B*F(1000001,1000000)
frac=lambda q:mp.mpf(F(q).numerator)/F(q).denominator
ps=[frac(x) for x in r['physical_polynomial']];MM=frac(M)
us=np.linspace(6,7,4001);hs=np.exp(-17/8*np.pi*(us-6));vs=[]
for u in us:
 X=mp.pi*float(u);h=mp.exp(-mp.mpf(17)/8*(X-6*mp.pi));v=mp.polyval(list(reversed(ps)),X)*mp.exp(-X)/MM/h;vs.append(float(v))
vs=np.array(vs)
def cost(lam):
 sig=np.max((vs+lam*hs)/2)
 return sig*sig/lam
fit=minimize_scalar(cost,bounds=(.01,100),method='bounded')
lam=fit.x;sig=max((vs+lam*hs)/2)
print({'lambda':lam,'sigma':sig,'cost_ratio':cost(lam),'argmax_u':float(us[np.argmax(vs+lam*hs)]),'at_one':cost(1),'max_f_at_6':float(vs[0])},flush=True)
(D/'origin_amplitude_screen.json').write_text(json.dumps({'lambda':lam,'sigma':sig,'cost_ratio':cost(lam),'at_one':cost(1),'certified':False},indent=2))
