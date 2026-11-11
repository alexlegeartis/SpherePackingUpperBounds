"""Generate/check the pinned comparison using standard-library Decimal arithmetic."""
from decimal import Decimal as D, localcontext
from fractions import Fraction
from pathlib import Path
import json, math, sys
R=Path(__file__).resolve().parent

def calculate(precision):
 with localcontext() as ctx:
  ctx.prec=precision
  def atan(x):
   return sum(((-1)**k*x**(2*k+1)/D(2*k+1) for k in range(240)),D(0))
  pi=16*atan(D(1)/5)-4*atan(D(1)/239)
  lows={6:pi**3/(48*D(3).sqrt()),45:pi**22*D(23)**22*D(23).sqrt()*D(math.factorial(23))/(D(2)**45*D(2).sqrt()*D(math.factorial(46))),46:(13*pi)**23/(D(math.factorial(23))*D(3)**46*D(3).sqrt()),47:pi**23*D(47)**23*D(47).sqrt()*D(math.factorial(24))/(D(2)**70*D(math.factorial(48))),17:pi**8*D(4)**9*D(math.factorial(9))/D(math.factorial(18))/16,18:pi**9/D(math.factorial(9))*(D(3)/4)**9,5:pi*pi/(15*D(2).sqrt()),4:pi*pi/16,21:pi**10*D(4)**11*D(math.factorial(11))/D(math.factorial(22))*D(43)**10*D(43).sqrt()/(D(2)**41*D(3)**11*D(3).sqrt()),19:pi**9*D(4)**10*D(math.factorial(10))/D(math.factorial(20))*D(13)**9*D(13).sqrt()/(D(3)**9*D(5)**10*D(5).sqrt()),20:pi**10/D(math.factorial(10))*D(3)**20/(D(2)**10*D(5)**10*D(5).sqrt()),48:(3*pi/2)**24/D(math.factorial(24))}
  rows=json.loads((R/'RESULTS.json').read_text())
  for row in rows:
   n=row['dimension']; L=lows[n]; U=D({6:'.4102905',45:'.0000013300905665',46:'.0000009295151556',47:'.0000006490757338',17:'.0184640903350649',18:'.0134853404450862',5:'.5122645',4:'.63610277',21:'.0051596603948176',19:'.0098179551395438',20:'.0071270536033763',48:'.0000004529067791'}[n]); N=D(row['strict_density_upper_bound'])
   assert 0<L<N<U
   row.update(known_lower_density=format(L,'.20g'),previous_upper_density=str(U),upper_decrease_percent=format(100*(U-N)/U,'.6f'),gap_decrease_percent=format(100*(U-N)/(U-L),'.6f'),old_gap=format(U-L,'.20g'),new_gap=format(N-L,'.20g'))
   if n==6:
    row.update(previous_upper_status='Paper-reported value; no independent replay',previous_upper_source='https://arxiv.org/html/2512.04829v2#S3.T1',external_neighbor_count_source='https://arxiv.org/html/1206.2608#S6.T2')
   if n==4:
    row.update(previous_upper_status='Paper-reported value printed to 8 decimal places; no independent replay',previous_upper_source='https://arxiv.org/html/2512.04829v2#S3.T1',previous_upper_reviewed='2026-11-10',historical_upper_density='0.6361073321551329',historical_upper_source='https://arxiv.org/abs/2206.15373')
   if n==5:
    row.update(previous_upper_status='Paper-reported value printed to 7 decimal places; no independent replay',previous_upper_source='https://arxiv.org/html/2512.04829v2#S3.T1',previous_upper_reviewed='2026-11-10',historical_upper_density='0.5126451306253027',historical_upper_source='https://arxiv.org/abs/2206.15373',external_neighbor_count_source='https://arxiv.org/html/1206.2608#S6.T2')
   # Independently check upper percentage with rational arithmetic.
   exact=100*(Fraction(U)-Fraction(N))/Fraction(U)
   assert abs(D(exact.numerator)/D(exact.denominator)-100*(U-N)/U)<D(10)**(-precision+5)
   assert abs((U-L)-(N-L)-(U-N))<D(10)**(-precision+5)
  return rows
rows=calculate(70)
assert rows==calculate(100), 'Displayed digits must be precision-stable'
lower_sources={6:'https://cohn.mit.edu/sphere-packing/',45:'https://arxiv.org/abs/2505.02394v1',46:'https://doi.org/10.1007/BF01232378',47:'https://arxiv.org/abs/2505.02394v1',17:'https://doi.org/10.4153/CJM-1964-065-1',18:'https://doi.org/10.1023/A:1008733715001',5:'https://cohn.mit.edu/sphere-packing/',4:'https://cohn.mit.edu/sphere-packing/',21:'https://arxiv.org/abs/2505.02394v1',19:'https://arxiv.org/abs/2505.02394v1',20:'https://arxiv.org/abs/2505.02394v1',48:'https://doi.org/10.4153/CJM-1971-081-3'}
upper_source={6:'https://arxiv.org/html/2512.04829v2#S3.T1',45:'https://arxiv.org/abs/2006.02560v2',46:'https://arxiv.org/abs/2006.02560v2',47:'https://arxiv.org/abs/2006.02560v2',17:'https://arxiv.org/abs/2006.02560v2',18:'https://arxiv.org/abs/2006.02560v2',5:'https://arxiv.org/html/2512.04829v2#S3.T1',4:'https://arxiv.org/html/2512.04829v2#S3.T1',19:'https://arxiv.org/abs/2006.02560v2',20:'https://arxiv.org/abs/2006.02560v2',21:'https://arxiv.org/abs/2006.02560v2',48:'https://arxiv.org/abs/2006.02560v2'}
header='| Dimension | Known lower density L | Older upper U | New upper N | Upper decrease | Additive gap decrease | Proof and details |\n| --- | --- | --- | --- | --- | --- | --- |\n'
def display_density(value):
 return '≈' + format(D(value), '.5g')
def display_percent(value):
 return format(D(value), '.2f') + '%'
assert display_density('0.006979502') == '≈0.0069795'
assert display_percent('4.755298') == '4.76%'
table=header+''.join(f"| {x['dimension']} | [{display_density(x['known_lower_density'])}]({lower_sources[x['dimension']]}) | [{display_density(x['previous_upper_density'])}]({upper_source[x['dimension']]}) | {display_density(x['strict_density_upper_bound'])} | {display_percent(x['upper_decrease_percent'])} | {display_percent(x['gap_decrease_percent'])} | [{x['dimension']}D proof]({x['package']}/) · [Sources & formulas](COMPARISONS.md) |\n" for x in rows)
display_note='Densities shown approximately to five significant digits; percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. The 4D, 5D and 6D baselines are pinned paper-reported decimal values, not independently replayed here. '
text='# Density bounds and pinned comparisons\n\n'+table+'\n'+display_note+'See [sources and conventions](COMPARISONS.md).\n'
outputs={'RESULTS.json':json.dumps(rows,indent=2)+'\n','RESULTS.md':text}
readme=(R/'README.md').read_text();a=readme.index('## Registered results');b=readme.index('These certificates use',a)
outputs['README.md']=readme[:a]+'## Registered results\n\n'+table+'\n'+display_note+'\n\n'+readme[b:]
for name,content in outputs.items():
 if '--write' in sys.argv:(R/name).write_text(content)
 else:assert (R/name).read_text()==content,name+' needs regeneration: python comparison.py --write'
print('PASS: ordering, independent rational percentages, gap identities, 70/100-digit agreement and generated tables')
