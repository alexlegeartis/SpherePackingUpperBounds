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
  lows={20:pi**10/D(math.factorial(10))*D(3)**20/(D(2)**10*D(5)**10*D(5).sqrt()),48:(3*pi/2)**24/D(math.factorial(24))}
  rows=json.loads((R/'RESULTS.json').read_text())
  for row in rows:
   n=row['dimension']; L=lows[n]; U=D({20:'.0071270536033763',48:'.0000004529067791'}[n]); N=D(row['strict_density_upper_bound'])
   assert 0<L<N<U
   row.update(known_lower_density=format(L,'.20g'),previous_upper_density=str(U),upper_decrease_percent=format(100*(U-N)/U,'.6f'),gap_decrease_percent=format(100*(U-N)/(U-L),'.6f'),old_gap=format(U-L,'.20g'),new_gap=format(N-L,'.20g'))
   # Independently check upper percentage with rational arithmetic.
   exact=100*(Fraction(U)-Fraction(N))/Fraction(U)
   assert abs(D(exact.numerator)/D(exact.denominator)-100*(U-N)/U)<D(10)**(-precision+5)
   assert abs((U-L)-(N-L)-(U-N))<D(10)**(-precision+5)
  return rows
rows=calculate(70)
assert rows==calculate(100), 'Displayed digits must be precision-stable'
header='| Dimension | Known lower density L | Older upper U | New strict upper N | Upper decrease | Additive gap decrease |\n| --- | --- | --- | --- | --- | --- |\n'
table=header+''.join(f"| {x['dimension']} | {x['known_lower_density']} | {x['previous_upper_density']} | {x['strict_density_upper_bound']} | {x['upper_decrease_percent']}% | {x['gap_decrease_percent']}% |\n" for x in rows)
text='# Density bounds and pinned comparisons\n\n'+table+'\nLower densities are approximate evaluations of exact formulas; U is the upward-rounded table entry and N the repository’s safely rounded strict bound. Percentages use N, not a tighter internal objective. See [sources and conventions](COMPARISONS.md).\n'
outputs={'RESULTS.json':json.dumps(rows,indent=2)+'\n','RESULTS.md':text}
readme=(R/'README.md').read_text();a=readme.index('## Registered results');b=readme.index('Both use',a)
outputs['README.md']=readme[:a]+'## Registered results\n\n'+table+'\n[References, exact lower-bound formulas and comparison conventions](COMPARISONS.md). [20D proof](verifications/dim20/) · [48D proof](verifications/dim48/).\n\n'+readme[b:]
for name,content in outputs.items():
 if '--write' in sys.argv:(R/name).write_text(content)
 else:assert (R/name).read_text()==content,name+' needs regeneration: python comparison.py --write'
print('PASS: ordering, independent rational percentages, gap identities, 70/100-digit agreement and generated tables')
