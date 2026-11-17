from flint import arb,ctx
from pathlib import Path
import json
ctx.prec=256
out={}
for n,det in [(83,4),(84,5),(85,6)]:
 d=arb(2)**(arb(n)/2)/arb(det).sqrt()*arb.pi()**(arb(n)/2)/(arb(n)/2+1).gamma()
 out[str(n)]={'density_interval':str(d),'center_density_squared':str(2**n//det)if 2**n%det==0 else f'{2**n}/{det}'}
print(out);Path(__file__).with_name('densities.json').write_text(json.dumps(out,indent=2))
