import json,math,sys
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq,fmpq_poly
ROOT=Path(__file__).resolve().parent;D=ROOT/'inputs';n=int(sys.argv[1]);d=json.loads((D/f'angular{n}_a40_scaled_candidate.json').read_text());assert d['dimension']==n
cs=[F(v).limit_denominator(10**12) for v in d['coefficients']];cs[0]-=F(1,1000);cs=[x/cs[0] for x in cs];assert min(cs)>=0 and cs[0]==1
ff=lambda x:fmpq(x.numerator,x.denominator);fr=lambda x:F(int(x.numerator),int(x.denominator));t=fmpq_poly([0,1]);G=[fmpq_poly([1]),t]
for l in range(1,len(cs)-1):G.append(((2*l+n-2)*t*G[-1]-l*G[-2])/(l+n-2))
Q=sum((ff(c)*g for c,g in zip(cs,G)),fmpq_poly());B=sum(cs);out={'dimension':n,'degree':Q.degree(),'alpha':40,'cutoff':'5001/10000','gegenbauer_coefficients':[str(v) for v in cs],'polynomial_coefficients':[str(fr(Q[i])) for i in range(20)],'B_exact':str(B),'B_float':float(B),'subtracted_epsilon':'1/1000','negative_cap_delta':'1/100000'}
S=ROOT/'work';S.mkdir(exist_ok=True);(S/'cyclic_angular40_exact.json').write_text(json.dumps(out,indent=2));print('B',float(B))
