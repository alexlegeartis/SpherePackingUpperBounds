from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import json,gzip,sys
sys.set_int_max_str_digits(0)
R=Path(__file__).resolve().parents[2];rows=json.loads((R/'LATTICE_RESULTS.json').read_text());seen=set();report=[]
for v in rows:
 n=v['dimension'];assert n not in seen;seen.add(n);p=R/v['precise_bound_evidence'];z=json.loads(gzip.open(p,'rt').read()if p.suffix=='.gz'else p.read_text());bound=F(v['certified_density_upper_bound'])
 if n<=35:assert F(z.get('safe_density',z['density']))<=bound
 elif n<=45:
  q=next(t for t in z if t['dimension']==n);assert F(q.get('upper',q.get('verified_upper')))==bound
 elif n==46:assert F(z['safe_bound'])==bound
 else:assert F(z.get('safe_upper',z.get('safe_bound')))==bound
 assert 'Kravatskiy'not in v['previous_method_source'];external=F(v['previous_lattice_upper']);lower=F(v['named_lattice_lower_comparison']);assert 0<lower<bound<external
 with localcontext()as ctx:
  ctx.prec=60
  dec=lambda x:D(x.numerator)/D(x.denominator)
  gain=dec(100*(external-bound)/external);gap=dec(100*(external-bound)/(external-lower));assert abs(gain-D(v['upper_decrease_percent']))<D('1e-48');assert abs(gap-D(v['gap_decrease_percent']))<D('1e-48')
 if n in [41,42]:assert v['external_comparator_certified']is False and v['previous_source_certificate_upper']is None and v['external_numerical_estimate']==v['previous_lattice_upper']
 else:assert v['external_comparator_certified']is True and F(v['previous_source_certificate_upper'])==external
 report.append({'dimension':n,'upper':str(bound),'evidence':v['precise_bound_evidence'],'comparison_certified':v['external_comparator_certified']})
assert len(rows)==34
print('PASS34 precise published bound references, exact comparison arithmetic and certification labels')
(Path(__file__).parent/'result.json').write_text(json.dumps(report,indent=2)+'\n')
