from pathlib import Path
import json
from flint import fmpq_poly,fmpq,ctx
ctx.prec=1024;D=Path(__file__).parent;z=json.load(open(D/'repaired.json'));P=fmpq_poly(list(map(fmpq,z['P'])));roots=P.complex_roots();assert sum(m for r,m in roots)==P.degree()
assert all((not r.imag.contains(0))or(r.real<180) for r,m in roots);assert P(180)<0
out=dict(status='PASS',all_roots_accounted=P.degree(),negative_from=180,leading_negative=P[P.degree()]<0);(D/'tail_result.json').write_text(json.dumps(out));print(out)
