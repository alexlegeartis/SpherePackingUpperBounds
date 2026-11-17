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
  lows={9:pi**4*D(2).sqrt()/D(945),10:pi**5/D(3072),16:pi**8/D(645120),11:D(144)*pi**5*D(math.factorial(6))/D(math.factorial(12)),12:pi**6/(D(27)*D(math.factorial(6))),13:D(576)*pi**6*D(math.factorial(7))/D(math.factorial(14)),14:pi**7/(D(16)*D(3).sqrt()*D(math.factorial(7))),15:D(2)**12*pi**7*D(math.factorial(8))/(D(2).sqrt()*D(math.factorial(16))),37:D(2)**40*D(2).sqrt()*pi**18*D(math.factorial(19))/D(math.factorial(38)),35:D(2)**37*D(2).sqrt()*pi**17*D(math.factorial(18))/D(math.factorial(36)),33:D(2)**9*pi**16*D(3)**16*D(3).sqrt()*D(math.factorial(17))/D(math.factorial(34)),31:D(2)**8*D(2).sqrt()*pi**15*D(3)**15*D(math.factorial(16))/D(math.factorial(32)),30:pi**15/D(math.factorial(15)),29:D(2)**30*pi**14*D(math.factorial(15))/(D(math.factorial(30))*D(2).sqrt()),23:pi**11*D(math.factorial(12))*D(23)**11*D(23).sqrt()/(D(2)**10*D(math.factorial(24))*D(3)**12),22:pi**11*D(11)**11/(D(math.factorial(11))*D(2)**23*D(3)**10*D(3).sqrt()),27:D(2)**28*pi**13*D(math.factorial(14))/(D(math.factorial(28))*D(2).sqrt()),26:pi**13/(D(math.factorial(13))*D(3).sqrt()),25:D(2)**26*pi**12*D(math.factorial(13))/(D(math.factorial(26))*D(2).sqrt()),28:pi**14/D(math.factorial(14)),32:pi**16*D(3)**16/(D(math.factorial(16))*D(2)**24),34:pi**17*D(3)**16*D(3).sqrt()/(D(math.factorial(17))*D(2)**25),36:D(16777216)*pi**18/(D(3176523)*D(math.factorial(18))),38:pi**19*D(3)**15*D(3).sqrt()/(D(math.factorial(19))*D(2)**19*5),39:D(2)**20*pi**19*D(3)**16*D(math.factorial(20))/(D(5).sqrt()*D(math.factorial(40))),40:pi**20*D(3)**16/(D(2)**20*D(math.factorial(20))),41:D(2)**21*pi**20*D(3)**17*D(math.factorial(21))/D(math.factorial(42)),42:pi**21*D(3)**17*D(3).sqrt()/(D(2)**21*D(math.factorial(21))),43:pi**21*D(137)**21*D(137).sqrt()*D(math.factorial(22))/(D(2)**67*D(2).sqrt()*D(3)**21*D(math.factorial(44))),7:pi**3/105,44:(157*pi)**22/(D(math.factorial(22))*D(2)**22*D(5)**21*D(5).sqrt()*D(11)**23),6:pi**3/(48*D(3).sqrt()),45:pi**22*D(23)**22*D(23).sqrt()*D(math.factorial(23))/(D(2)**45*D(2).sqrt()*D(math.factorial(46))),46:(13*pi)**23/(D(math.factorial(23))*D(3)**46*D(3).sqrt()),47:pi**23*D(47)**23*D(47).sqrt()*D(math.factorial(24))/(D(2)**70*D(math.factorial(48))),17:pi**8*D(4)**9*D(math.factorial(9))/D(math.factorial(18))/16,18:pi**9/D(math.factorial(9))*(D(3)/4)**9,5:pi*pi/(15*D(2).sqrt()),4:pi*pi/16,21:pi**10*D(4)**11*D(math.factorial(11))/D(math.factorial(22))*D(43)**10*D(43).sqrt()/(D(2)**41*D(3)**11*D(3).sqrt()),19:pi**9*D(4)**10*D(math.factorial(10))/D(math.factorial(20))*D(13)**9*D(13).sqrt()/(D(3)**9*D(5)**10*D(5).sqrt()),20:pi**10/D(math.factorial(10))*D(3)**20/(D(2)**10*D(5)**10*D(5).sqrt()),48:(3*pi/2)**24/D(math.factorial(24))}
  rows=json.loads((R/'RESULTS.json').read_text())
  for row in rows:
   n=row['dimension']; L=lows[n]; U=D({9:'0.19099063',10:'0.14340271',16:'0.02492121',11:'0.10672201',12:'0.0797066',13:'0.06012103',14:'0.04501031',15:'0.03372585',37:'2.26586900106e-05',35:'4.55743843107e-05',33:'9.12356039023e-05',31:'0.0001817083813917',30:'0.0002559028743732',29:'0.0003598581852089',23:'0.0026842798864291',22:'0.0037259419689206',27:'0.0007082297958617',26:'0.0009910238892216',25:'0.0013841907222857',28:'0.0005052542161057',32:'0.0001288432887595',34:'6.45221967438e-05',36:'0.0000321530553313',38:'1.59506499105e-05',39:'1.12168687009e-05',40:'0.0000078801051697',41:'.0000055306464395',42:'.0000038780970907',43:'.0000027169074727',7:'.32111475',44:'.0000019017703144',6:'.4102905',45:'.0000013300905665',46:'.0000009295151556',47:'.0000006490757338',17:'.0184640903350649',18:'.0134853404450862',5:'.5122645',4:'.63610277',21:'.0051596603948176',19:'.0098179551395438',20:'.0071270536033763',48:'.0000004529067791'}[n]); N=D(row['strict_density_upper_bound'])
   # Dimension 4 now compares against the fully replayed, pinned Felpix proof.
   if n==4:U=D('0.629521681653')
   assert 0<L<N<U
   row.update(known_lower_density=format(L,'.20g'),previous_upper_density=str(U),upper_decrease_percent=format(100*(U-N)/U,'.6f'),gap_decrease_percent=format(100*(U-N)/(U-L),'.6f'),old_gap=format(U-L,'.20g'),new_gap=format(N-L,'.20g'))
   if n==6:
    row.update(previous_upper_status='Paper-reported value; no independent replay',previous_upper_source='https://arxiv.org/html/2512.04829v2#S3.T1',external_neighbor_count_source='https://arxiv.org/html/1206.2608#S6.T2')
   if n==4:
    row.update(previous_upper_status='Pinned Felpix-Studios global certificate; complete full-mode replay including all 112 cyclic branches',previous_upper_source='https://github.com/Felpix-Studios/sphere-packing/tree/bb97550ba47ad92bacfee567843888be16895a23/preprints/D4',previous_upper_reviewed='2026-11-15',historical_upper_density='0.6361073321551329',historical_upper_source='https://arxiv.org/abs/2206.15373')
   if n==5:
    row.update(previous_upper_status='Paper-reported value printed to 7 decimal places; no independent replay',previous_upper_source='https://arxiv.org/html/2512.04829v2#S3.T1',previous_upper_reviewed='2026-11-10',historical_upper_density='0.5126451306253027',historical_upper_source='https://arxiv.org/abs/2206.15373',external_neighbor_count_source='https://arxiv.org/html/1206.2608#S6.T2')
   # Independently check upper percentage with rational arithmetic.
   exact=100*(Fraction(U)-Fraction(N))/Fraction(U)
   assert abs(D(exact.numerator)/D(exact.denominator)-100*(U-N)/U)<D(10)**(-precision+5)
   assert abs((U-L)-(N-L)-(U-N))<D(10)**(-precision+5)
  return rows
rows=calculate(70)
assert rows==calculate(100), 'Displayed digits must be precision-stable'
lower_sources={9:'https://cohn.mit.edu/sphere-packing/',10:'https://cohn.mit.edu/sphere-packing/',16:'https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/BW16.html',11:'https://cohn.mit.edu/sphere-packing/',12:'https://cohn.mit.edu/sphere-packing/',13:'https://cohn.mit.edu/sphere-packing/',14:'https://cohn.mit.edu/sphere-packing/',15:'https://cohn.mit.edu/sphere-packing/',37:'https://cohn.mit.edu/sphere-packing/',35:'https://cohn.mit.edu/sphere-packing/',33:'https://cohn.mit.edu/sphere-packing/',31:'https://cohn.mit.edu/sphere-packing/',30:'https://cohn.mit.edu/sphere-packing/',29:'https://cohn.mit.edu/sphere-packing/',23:'https://cohn.mit.edu/sphere-packing/',22:'https://cohn.mit.edu/sphere-packing/',27:'https://cohn.mit.edu/sphere-packing/',26:'https://cohn.mit.edu/sphere-packing/',25:'https://cohn.mit.edu/sphere-packing/',28:'https://cohn.mit.edu/sphere-packing/',32:'https://cohn.mit.edu/sphere-packing/',34:'https://cohn.mit.edu/sphere-packing/',36:'constructions/dim36/',38:'https://arxiv.org/html/2607.20359v4',39:'https://arxiv.org/html/2607.20359v4',40:'https://arxiv.org/html/2607.20359v4',41:'https://arxiv.org/html/2607.20359v4',42:'https://arxiv.org/html/2607.20359v4',43:'https://arxiv.org/html/2607.20359v4',7:'https://cohn.mit.edu/sphere-packing/',44:'https://arxiv.org/abs/2505.02394v1',6:'https://cohn.mit.edu/sphere-packing/',45:'https://arxiv.org/abs/2505.02394v1',46:'https://doi.org/10.1007/BF01232378',47:'https://arxiv.org/abs/2505.02394v1',17:'https://doi.org/10.4153/CJM-1964-065-1',18:'https://doi.org/10.1023/A:1008733715001',5:'https://cohn.mit.edu/sphere-packing/',4:'https://cohn.mit.edu/sphere-packing/',21:'https://arxiv.org/abs/2505.02394v1',19:'https://arxiv.org/abs/2505.02394v1',20:'https://arxiv.org/abs/2505.02394v1',48:'https://doi.org/10.4153/CJM-1971-081-3'}
upper_source={9:'https://arxiv.org/html/2512.04829v2#S3.T1',10:'https://arxiv.org/html/2512.04829v2#S3.T1',16:'https://arxiv.org/html/2512.04829v2#S3.T1',11:'https://arxiv.org/html/2512.04829v2#S3.T1',12:'https://arxiv.org/html/2512.04829v2#S3.T1',13:'https://arxiv.org/html/2512.04829v2#S3.T1',14:'https://arxiv.org/html/2512.04829v2#S3.T1',15:'https://arxiv.org/html/2512.04829v2#S3.T1',37:"https://arxiv.org/abs/2006.02560v2",35:"https://arxiv.org/abs/2006.02560v2",33:"https://arxiv.org/abs/2006.02560v2",31:"https://arxiv.org/abs/2006.02560v2",30:"https://arxiv.org/abs/2006.02560v2",29:"https://arxiv.org/abs/2006.02560v2",23:"https://arxiv.org/abs/2006.02560v2",22:"https://arxiv.org/abs/2006.02560v2",27:"https://arxiv.org/abs/2006.02560v2",26:"https://arxiv.org/abs/2006.02560v2",25:"https://arxiv.org/abs/2006.02560v2",28:"https://arxiv.org/abs/2006.02560v2",32:"https://arxiv.org/abs/2006.02560v2",34:"https://arxiv.org/abs/2006.02560v2",36:"https://arxiv.org/abs/2006.02560v2",38:"https://arxiv.org/abs/2006.02560v2",39:"https://arxiv.org/abs/2006.02560v2",40:"https://arxiv.org/abs/2006.02560v2",41:'https://arxiv.org/abs/2006.02560v2',42:'https://arxiv.org/abs/2006.02560v2',43:'https://arxiv.org/abs/2006.02560v2',7:'https://arxiv.org/html/2512.04829v2#S3.T1',44:'https://arxiv.org/abs/2006.02560v2',6:'https://arxiv.org/html/2512.04829v2#S3.T1',45:'https://arxiv.org/abs/2006.02560v2',46:'https://arxiv.org/abs/2006.02560v2',47:'https://arxiv.org/abs/2006.02560v2',17:'https://arxiv.org/abs/2006.02560v2',18:'https://arxiv.org/abs/2006.02560v2',5:'https://arxiv.org/html/2512.04829v2#S3.T1',4:'https://arxiv.org/html/2512.04829v2#S3.T1',19:'https://arxiv.org/abs/2006.02560v2',20:'https://arxiv.org/abs/2006.02560v2',21:'https://arxiv.org/abs/2006.02560v2',48:'https://arxiv.org/abs/2006.02560v2'}
upper_source[4]='https://github.com/Felpix-Studios/sphere-packing/tree/bb97550ba47ad92bacfee567843888be16895a23/preprints/D4'
# The method/source beside U describes the recorded comparison, not our method.
# Keep this attribution in the registry as well as the generated Markdown.
comparison_families = {
 'felpix': {
  'method': 'Scalar/kernel + count and far-pair refinements',
  'label': 'Felpix-Studios, pinned bb97550',
  'title': 'An improved upper bound for sphere packing density in four dimensions',
  'locator': 'preprints/D4 at commit bb97550ba47ad92bacfee567843888be16895a23',
  'kind': 'pinned_repository_certificate',
 },
 'three_point': {
  'method': 'Three-point SDP',
  'label': 'Tutunov et al. (2025), Table 1',
  'title': 'Model-Based and Sample-Efficient AI-Assisted Math Discovery in Sphere Packing',
  'locator': 'arXiv:2512.04829v2, Table 1',
  'kind': 'paper_reported_decimal',
 },
 'linear_programming': {
  'method': 'Cohn–Elkies LP / modular bootstrap',
  'label': 'Afkhami-Jeddi et al. (2020)',
  'title': 'High-dimensional sphere packing and the modular bootstrap',
  'locator': 'arXiv:2006.02560v2; numerical comparator retained from the pinned Cohn table',
  'kind': 'published_table_decimal',
 },
}
for row in rows:
 n=row['dimension']
 family='felpix' if n==4 else 'three_point' if n<=16 else 'linear_programming'
 source=comparison_families[family]
 assert row.get('previous_upper_source',upper_source[n])==upper_source[n]
 row.update(previous_upper_method=source['method'],
            previous_upper_source=upper_source[n],
            previous_upper_source_label=source['label'],
            previous_upper_source_title=source['title'],
            previous_upper_source_locator=source['locator'],
            previous_upper_source_kind=source['kind'])
# Central mathematical mechanisms, not just verification-package labels.
# Shared labels deliberately identify the same method across dimensions.
proof_ideas = {9:'Signed-hinge cyclic profile + average-neighbor count',10:'Log-hinge cyclic profile + root-cut radial repair',16:'Widened shell count + coupled log-hinge radial kernel',11:'Quartic radial damping + all-dilation cyclic reduction',12:'Quartic radial damping + all-dilation cyclic reduction',13:'Quartic radial damping + all-dilation cyclic reduction',14:'Quartic radial damping + all-dilation cyclic reduction',15:'Quartic radial damping + all-dilation cyclic reduction',37:"Dimension-specific angular kernel + Gaussian radial optimization",35:"Dimension-specific angular kernel + Gaussian radial optimization",33:"Dimension-specific angular kernel + Gaussian radial optimization",31:"Dimension-specific angular kernel + Gaussian radial optimization",30:"Dimension-specific angular kernel + Gaussian radial optimization",29:"Dimension-specific angular kernel + Gaussian radial optimization",23:"Dimension-specific angular kernel + Gaussian radial optimization",22:"Dimension-specific angular kernel + Gaussian radial optimization",27:"Dimension-specific angular kernel + Gaussian radial optimization",26:"Dimension-specific angular kernel + Gaussian radial optimization",
 4: 'Degree16 far-pair polynomial + unused angular Gram credits',
 5: 'Spherical lens caps + long-edge PSD accounting',
 6: 'Neighbor-count cuts + Schur-saturated cyclic kernel',
 7: 'Neighbor-count cuts + Schur-saturated cyclic kernel',
 17: 'Mixed-Gaussian tail repair + cyclic origin saturation',
 18: 'Solid-harmonic cyclic kernel + Schur origin saturation',
 19: 'Solid-harmonic cyclic kernel + Schur origin saturation',
 20: 'Solid-harmonic cyclic kernel + Schur origin saturation',
 21: 'Dimension-specific angular kernel + Gaussian radial optimization',
 25: 'Dimension-specific angular kernel + Gaussian radial optimization',
 28: 'Dimension-specific angular kernel + Gaussian radial optimization',
 32: 'Dimension-specific angular kernel + Gaussian radial optimization',
 34: 'Dimension-specific angular kernel + Gaussian radial optimization',
 36: 'Dimension-specific angular kernel + Gaussian radial optimization',
 38: 'Dimension-specific angular kernel + Gaussian radial optimization',
 39: 'Dimension-specific angular kernel + Gaussian radial optimization',
 40: 'Dimension-specific angular kernel + Gaussian radial optimization',
 41: 'Dimension-specific angular kernel + Gaussian radial optimization',
 42: 'Dimension-specific angular kernel + Gaussian radial optimization',
 43: 'Dimension-specific angular kernel + Gaussian radial optimization',
 44: 'Dimension-specific angular kernel + Gaussian radial optimization',
 45: 'Extended pair coupling + degree47 radial correction',
 46: 'Extended pair coupling + degree63 radial correction',
 47: 'Extended pair coupling + degree47 radial correction',
 48: 'Gaussian cyclic kernel + coupled radial optimization',
}
for n in [22, 25, 32, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]:proof_ideas[n]='Signed logarithmic profile + degree95 radial enrichment'
for n in [26, 27, 28, 29, 30, 31, 33, 34]:proof_ideas[n]='Signed logarithmic profile + degree95 radial enrichment'
for n in [23]:proof_ideas[n]='Signed logarithmic profile + degree95 radial enrichment'
for n in [17,18,19,20,21,48]:proof_ideas[n]='Signed logarithmic profile + degree95 radial enrichment'
for n in [45,48]:proof_ideas[n]='Reoptimized nonnegative Gegenbauer cyclic kernel + signed logarithmic profile + retained Fourier witness'
assert set(proof_ideas) == {row["dimension"] for row in rows}
header='| Dimension | Known lower density L | Recorded upper U | Comparison method / source | Our upper N | Upper decrease | Additive gap decrease | Our method / proof |\n| --- | --- | --- | --- | --- | --- | --- | --- |\n'
def display_density(value):
 return '≈' + format(D(value), '.5g')
def display_percent(value,dimension=None):
 return format(D(value), '.6f' if dimension==4 else '.2f') + '%'
def display_upper(row,field):
 if row['dimension']==4:
  return ('<' if field=='strict_density_upper_bound' else '')+row[field]
 return display_density(row[field])
assert display_density('0.006979502') == '≈0.0069795'
assert display_percent('4.755298') == '4.76%'
# Only verified lower-bound constructions originating in this repository.
this_work_lower_dimensions = {36}
def display_lower(row):
 label = display_density(row['known_lower_density'])
 if row['dimension'] in this_work_lower_dimensions:
  label += ' (this work)'
 return label
assert display_lower({'dimension': 36, 'known_lower_density': '7.33034778317e-7'}) == '≈7.3303e-7 (this work)'
assert display_lower({'dimension': 16, 'known_lower_density': '0.014708'}) == '≈0.014708'
table=header+''.join(f"| {x['dimension']} | [{display_lower(x)}]({lower_sources[x['dimension']]}) | {display_upper(x,'previous_upper_density')} | [{x['previous_upper_method']} — {x['previous_upper_source_label']}]({x['previous_upper_source']}) | {display_upper(x,'strict_density_upper_bound')} | {display_percent(x['upper_decrease_percent'],x['dimension'])} | {display_percent(x['gap_decrease_percent'],x['dimension'])} | [{proof_ideas[x['dimension']]}]({x['package']}/) · [Sources & formulas](COMPARISONS.md) |\n" for x in rows)
comparison_note='**Recorded upper U** is the historical comparison used for the percentages in that row. **Comparison method / source** names the method behind U and links to its paper or pinned certificate; **Our method / proof** describes the new bound N. A recorded comparison is not a claim about the latest worldwide record. [Full source titles, versions and comparison scope](COMPARISONS.md#table-comparison-methods-and-sources).'
display_note='The 4D comparison uses the pinned Felpix-Studios certificate and its complete full-mode replay; extra digits show the small strict improvement. The 16D strict bound is 0.024900975, compared with the paper’s displayed 0.02492121. The full widened-shell count certificate is included and independently replayed. No worldwide-priority claim is made. Other densities are shown approximately to five significant digits and percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. The 5D, 6D, 7D and 9–16D baselines are pinned paper-reported decimal values, not independently replayed here.'
text='# Density bounds and pinned comparisons\n\n'+comparison_note+'\n\n'+table+'\n'+display_note+' See [sources and conventions](COMPARISONS.md).\n'
outputs={'RESULTS.json':json.dumps(rows,indent=2)+'\n','RESULTS.md':text}
readme=(R/'README.md').read_text();a=readme.index('## Registered results');b=readme.index('These certificates use',a)
compact_header='| Dimension / lower L | Recorded upper U / source | Our upper N / method | Upper decrease | Gap decrease |\n|---|---|---|---:|---:|\n'
compact_table=compact_header+''.join(f"| {x['dimension']}<br>**L:** [{display_lower(x)}]({lower_sources[x['dimension']]}) | {display_upper(x,'previous_upper_density')}<br>[{x['previous_upper_method']} — {x['previous_upper_source_label']}]({x['previous_upper_source']}) | {display_upper(x,'strict_density_upper_bound')}<br>[{proof_ideas[x['dimension']]}]({x['package']}/) · [Sources & formulas](COMPARISONS.md) | {display_percent(x['upper_decrease_percent'],x['dimension'])} | {display_percent(x['gap_decrease_percent'],x['dimension'])} |\n" for x in rows)
outputs['README.md']=readme[:a]+'## Registered results\n\n'+comparison_note+'\n\n'+compact_table+'\n'+display_note+'\n\n'+readme[b:]
for name,content in outputs.items():
 if '--write' in sys.argv:(R/name).write_text(content)
 else:assert (R/name).read_text()==content,name+' needs regeneration: python comparison.py --write'
print('PASS: ordering, independent rational percentages, gap identities, 70/100-digit agreement and generated tables')
