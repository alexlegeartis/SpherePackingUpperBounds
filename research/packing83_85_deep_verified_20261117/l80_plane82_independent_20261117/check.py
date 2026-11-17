from pathlib import Path
from fractions import Fraction as F
from math import factorial
import numpy as np,json,hashlib,decimal
from flint import fmpq_mat as Q
D=Path(__file__).parent;L=D.parent/'bachoc_nebe80_reconstruction_20261117';z=np.load(L/'lattice80.npz');M=Q(z['twice_root_coordinate_basis'].tolist());G=Q(z['gram'].tolist());Mi=M.inv();iso=json.loads((D/'isometries.json').read_text())['isometries'];items=json.loads((D/'deep_plane.json').read_text());classes=[]
assert len(items)==3
for it in items:
 case=it['source_case'];assert case in [101,108];p=D.parent/'l80_code_decoder_20261117'/f'case{case}.json';src=json.loads(p.read_text());assert Q(src['gram'])==G
 # Bind inherited source data to the already published independent proof inputs.
 expected={101: 'bced37fa26cf53923af5f5fc59379203939f3b96f3692e8945110c6c2351022b', 108: 'b7e5df20efbe94603a0cafbd123e36e80d3f41af3b55ad68546260fd185c907c'};assert hashlib.sha256(p.read_bytes()).hexdigest()==expected[case]
 U=Q(iso[it['isometry_index']]);bp=it['block_permutation'];assert sorted(bp)==list(range(10));mapped=[]
 for row in M.tolist():
  blocks=[None]*10
  for j in range(10):blocks[bp[j]]=(Q([row[8*j:8*j+8]])*U).tolist()[0]
  mapped.append(sum(blocks,[]))
 A=Q(mapped)*Mi;assert all(x.denominator==1 for x in A.entries()) and abs(A.det())==1 and A*G*A.transpose()==G
 lift=Q([src['r']])*A;assert lift==Q([it['lift']]);cl=[int(x)%2 for x in lift.entries()];assert cl==it['class'] and any(cl);classes.append(cl)
assert len(set(map(tuple,classes)))==3 and all((a+b+c)%2==0 for a,b,c in zip(*classes))
R=Q(classes[:2])/2;cross=G*R.transpose();bottom=R*G*R.transpose()+Q([[2,-1],[-1,2]]);rows=[G.tolist()[i]+cross.tolist()[i]for i in range(80)]+[(R*G).tolist()[i]+bottom.tolist()[i]for i in range(2)];H=Q(rows);c=json.loads((D/'candidate82_plane.json').read_text());assert H==Q(c['gram']) and H.det()==3
# A2 quadratic form is m²+n²+(m−n)², a positive even integer for nonzero(m,n).
# Hence >=2; for even-even nonzero coefficients >=8. Every other parity class has horizontal depth>=6.
delta=F(2)**82/3;assert delta==F(c['center_density_squared']);floor=F('9.17033181889403e-18');assert floor==F(c['ordinary_density_floor'])
at=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0));pi=16*at(5,201)-4*at(239,200);assert floor**2<delta*pi**82/factorial(41)**2
ext=F(7,2)**82/(147*41**12);assert ext==F(c['external_MW80_A2_center_density_squared']);decimal.getcontext().prec=60;gain=(decimal.Decimal((delta/ext).numerator)/decimal.Decimal((delta/ext).denominator)).sqrt()*100-100
out={'status':'COMPLETE INDEPENDENT EXACT TRANSFER/GEOMETRY/DENSITY AUDIT PASS','dimension':82,'density_floor':str(floor),'minimum_squared':8,'determinant':3,'source_depths':'Previously independently certified cases101 and108, bound to published exact input bytes','external_gain_percent':str(gain)};(D/'result.json').write_text(json.dumps(out,indent=2));print(out)
