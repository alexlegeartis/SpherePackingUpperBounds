import sys;sys.set_int_max_str_digits(0)
import json,math
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P
from pathlib import Path
D=Path(__file__).parent;z=json.load(open(D/'input_candidate.json'));r=json.load(open(D/'repaired.json'));a=json.load(open(D/'annular_final_certificate.json'));p=P(list(map(Q,r['P'])));t=p(P([180,1]));assert all(x<0 for x in t);val=(F(str(p[0]))+F(z['s'])+F(a['budget']))*F(17,4)**18/math.factorial(18);safe=F('2.727705e-5');assert val<safe
out={'status':'PRODUCER SCALAR CLOSURE PASS; INDEPENDENT GEOMETRY LINK AND REPLAY PENDING','tail':'all coefficients P(180+t) strictly negative','density_exact':str(val),'safe':str(safe),'annuli':len(a['rows']),'nodes':a['nodes']};(D/'scalar_result.json').write_text(json.dumps(out,indent=2));print('PASS',float(val),'safe',float(safe),'budget',float(F(a['budget'])))
