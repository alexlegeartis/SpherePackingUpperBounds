from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
P=Path(__file__).resolve().parent;R=P.parent.parent
release=json.loads((R/'research/packing_release8590_20261117/RESULTS.json').read_text())['rows']
old={x['dimension']:x for x in json.loads((R/'research/mw92_sections_20261115/results.json').read_text())['rows']}
b={85:F(2**80),86:F(2**84,3),87:F(2**85),88:F(2**88)}
kl={x['dimension']:F(x['center_density_squared_lower']) for x in json.loads((P/'kl_report.json').read_text())['sections']}
rows=[]
for r in release:
 n=r['dimension']; raw=json.loads((R/'research/packing_release8590_20261117/packing_resume_20261117'/r['file']).read_text());new=F(raw['density_squared'])
 if n<89:lo=hi=b[n]
 else:
  lo=hi=kl[n]
 with localcontext() as ctx:
  ctx.prec=80;ratio=new/hi;pct=((Decimal(ratio.numerator)/Decimal(ratio.denominator)).sqrt()-1)*100; k=int(pct*100)
 assert (1+F(k,10000))**2<new/hi
 assert new/hi<(1+F(k+1,10000))**2
 rows.append({'dimension':n,'external_center_density_squared_lower':str(lo),'external_center_density_squared_upper':str(hi),'derived':n!=88,'certified_lower_guarantee':True,'gain_percent_strict_lower':f'{k//100}.{k%100:02d}','our_density_floor':r['rigorous_occupied_density_lower'],'preceding_repository_gain_percent':r['gain_percent_vs_published_approx']})
out={'scope':'Strongest external-derived comparison certified in the documented17 November follow-up. Not an exhaustive global-priority claim. Derived values are our calculations from external theorems, not numbers claimed printed by those authors.','rows':rows}
if __name__=='__main__':
 (P/'comparison.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
