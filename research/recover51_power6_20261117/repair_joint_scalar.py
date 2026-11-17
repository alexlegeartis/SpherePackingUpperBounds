"""Build exact normalized Fourier polynomial from a NUMERICAL joint candidate."""
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P
import json,sys,shutil
sys.set_int_max_str_digits(0);D=Path(__file__).parent;c=Path(sys.argv[1]);v=json.load(open(c));n=v['dimension'];assert n==51
base=D/'base/scalar';raw=json.load(open(base/'repaired.json'));poly=P(list(map(Q,raw['Q'])));x=P([0,1]);L=[P([1]),P([Q(n,2),-2])]
for k in range(2,72):L.append(((2*k+Q(n,2)-2-2*x)*L[-1]-(k+Q(n,2)-2)*L[-2])/k)
for j,d in enumerate(v['radial_delta']):poly+=Q(str(F(format(d,'.16g'))))*(-1)**j*L[j]/L[j](0)
assert poly.degree()==71 and poly(0)>0 and poly[71]>0
poly/=poly(0);out=D/(c.stem+'_scalar_probe');out.mkdir(exist_ok=True);(out/'input_roots.json').write_text(json.dumps({'Q':list(map(str,poly))},indent=2));shutil.copy(base/'repair.py',out/'repair.py');(out/'STATUS.md').write_text('Scalar repair probe only. Kernel geometry and PSD may fail. NOT A BOUND.\n');print(out)
