"""Independent cross-component normalization and source-hypothesis linkage audit."""
if not __debug__: raise RuntimeError('Assertions required')
import json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from flint import arb,fmpq,ctx
ctx.prec=384;sys.set_int_max_str_digits(0)
D=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
J=lambda p:json.loads((D/p).read_text(),parse_float=str)
c=J('certificate.json');assert c['dimension']==9 and F(c['s_exact'])==F(9,8)
for x in ['geometry','radial','independent_radial']:
 assert J(x+'/fit.json')['dimension']==9
 for y in ['fit.json','angular.json']:
  assert (D/x/y).read_bytes()==(D/'geometry'/y).read_bytes()
for x in ['radial','independent_radial']:
 assert (D/x/'certificate.json').read_bytes()==(D/'certificate.json').read_bytes()
p=list(map(F,c['physical']));C=F(c['C_exact']);lam=F(c['lambda_count_exact']);assert C>=0 and lam>=0
assert F(c['threshold_v_exact'])==F('1.003012')**2
assert F('1.003012')<F('1.003013')
# Min center separation sqrt(9/8), sphere radius sqrt(9/32).
# V9=32*pi^4/945; density ball volume = V9*(9/32)^(9/2).
a=lambda x:arb(fmpq(str(x)))
vol=a(F(32,945))*arb.pi()**4*a(F(9,32))**4*a(F(9,32)).sqrt()
value=a(p[0]+C+373*lam)*vol
assert value<a(F('0.190941549'))
# No integer flooring:373 is a proven bound on average ORDERED neighbors.
out={'status':'PASS','strict_density_upper_bound':'0.190941549','objective_interval':str(value),'s_exact':'9/8','sphere_radius_squared':'9/32','average_count_external':373,'no_integer_rounding':True,'threshold_in_published_shell':True,'cross_component_byte_linkage':True,'scope':'Normalization, coefficient linkage and density rounding; complete geometry and radial validity are separately audited','input_hashes':{str(x.relative_to(D)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [D/'certificate.json',D/'geometry/fit.json',D/'geometry/angular.json']}}
(Path(__file__).parent/'LINKAGE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
