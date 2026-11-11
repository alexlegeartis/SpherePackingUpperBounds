from pathlib import Path
from flint import arb,ctx
from fractions import Fraction
import math,json
ctx.prec=512;pi=arb.pi();L=arb(2)**12*pi**7*math.factorial(8)/(arb(2).sqrt()*math.factorial(16));U=arb('0.03372585');N=arb('0.033671603');a=100*(U-N)/U;g=100*(U-N)/(U-L)
assert a>arb('.160846')and a<arb('.160848');assert g>arb('.321591')and g<arb('.321593')
out={'pass':True,'lower_formula':'2^12*pi^7*8!/(sqrt(2)*16!)','lower_source':'https://cohn.mit.edu/sphere-packing/','lower_source_status':'Indexed four months earlier, center density2^(-9/2); construction not reconstructed','upper_source':'https://arxiv.org/html/2512.04829v2#S3.T1','upper_status':'Live Table1 checked2026-11-11, printed .03372585; source certificate not replayed','lower_interval':str(L),'upper_decrease_percent_interval':str(a),'gap_decrease_percent_interval':str(g),'precision_bits':512}
(Path(__file__).resolve().parent/'COMPARISON_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
