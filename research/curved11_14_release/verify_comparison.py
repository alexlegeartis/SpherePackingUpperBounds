from pathlib import Path
from flint import arb,ctx
import math,json
ctx.prec=512;p=arb.pi();lows={11:144*p**5*math.factorial(6)/math.factorial(12),12:p**6/(27*math.factorial(6)),13:576*p**6*math.factorial(7)/math.factorial(14),14:p**7/(16*arb(3).sqrt()*math.factorial(7))};rows=json.loads((Path(__file__).resolve().parents[2]/'RESULTS.json').read_text());out=[]
for n,L in lows.items():
 r=next(x for x in rows if x['dimension']==n);U=arb(r['previous_upper_density']);N=arb(r['strict_density_upper_bound']);a=100*(U-N)/U;g=100*(U-N)/(U-L);assert abs(a-arb(r['upper_decrease_percent']))<arb('.00000051');assert abs(g-arb(r['gap_decrease_percent']))<arb('.00000051');out.append({'dimension':n,'pass':True,'lower_interval':str(L),'upper_decrease_percent_interval':str(a),'gap_decrease_percent_interval':str(g),'precision_bits':512})
(Path(__file__).resolve().parent/'COMPARISON_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print([(r['dimension'],r['upper_decrease_percent'],r['gap_decrease_percent'])for r in rows if r['dimension']in lows])
