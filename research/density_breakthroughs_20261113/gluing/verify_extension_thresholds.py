"""Exact density thresholds for lattice extensions of the current 36D construction."""
from fractions import Fraction as F
from pathlib import Path
import json,math
P=Path(__file__).resolve().parent
delta=F(16777216,3176523);h2=delta**2/16;r2=8-h2
assert h2<2
# delta37^2=2 delta36^2/h^2; reference37 has squared center density32.
assert 2*delta**2/h2==32
# Order2 offsets have2s inL, so thepure vertical vector(0,2h)forces4h^2>=8.
assert delta**2<32
out={'dimension36_center_density':str(delta),'dimension37_reference_center_density_squared':'32','necessary_extension_height_squared_upper':str(h2),'necessary_nearest_layer_hole_squared_lower':str(r2),'height_squared_upper_decimal':float(h2),'hole_squared_lower_decimal':float(r2),'order_two_extension_center_density_upper':str(delta),'order_two_can_beat_reference37':False,'scope':'Exact restrictions for stacking copies of the fixed 36D lattice; not a general 37D obstruction.'}
for f in sorted(P.glob('torsion_p*_seed*.json')):
 j=json.loads(f.read_text())['best'];Q=[[int(t) for t in row.split()]for row in (P.parents[2]/'constructions/dim36/Q_reduced.txt').read_text().splitlines()];a=j['shift_numerators'];p=j['denominator'];z=j['nearest_witness'];v=[a[i]-p*z[i]for i in range(36)];q=sum(v[i]*Q[i][k]*v[k]for i in range(36)for k in range(36));assert q==j['distance_squared_num'];assert F(q,20*p*p)<r2
out['saved_best_torsion_candidates_have_exact_rejection_witnesses']=True
(P/'extension_thresholds.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
