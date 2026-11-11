"""Evaluate published packing formulas, not a new construction proof."""
from pathlib import Path
import json,math
from flint import arb,ctx
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=256;D=Path(__file__).resolve().parent
source='https://cohn.mit.edu/sphere-packing/'
upper={22:'0.0037259419689206',23:'0.0026842798864291'}
new={22:'0.003653383161',23:'0.002642540089'}
center={22:arb(11)**11/(arb(2)**23*arb(3)**10*arb(3).sqrt()),23:arb(23)**11*arb(23).sqrt()/(arb(2)**34*arb(3)**12)}
formulas={22:'2^(-23) * 3^(-21/2) * 11^11',23:'2^(-34) * 3^(-12) * 23^(23/2)'}
prefix={22:'0.0024510340441211',23:'0.0019071945388448'}
rows=[]
for n in [22,23]:
 volume=arb.pi()**(arb(n)/2)/(arb(n)/2+1).gamma();density=center[n]*volume
 assert density>arb(prefix[n]) and density<arb(prefix[n])+arb('0.0000000000000001')
 old=arb(upper[n]);current=arb(new[n]);assert density<current<old
 row={'dimension':n,'source_url':source,'retrieved':'2026-11-11','source_status':'Search-indexed Cohn table, crawled four months before retrieval; direct page returned502.',
 'center_density_formula':formulas[n],'packing_density_formula':f'pi^({n}/2)/Gamma({n}/2+1) * ({formulas[n]})',
 'center_density_interval':str(center[n]),'known_packing_density_interval':str(density),
 'displayed_density_prefix':prefix[n],'periodic_centers_per_cell':3 if n==22 else 2,
 'table_reference_indices':'[10,1]'if n==22 else'[6,1]',
 'pinned_old_upper':upper[n],'new_strict_upper':new[n],
 'upper_reduction_percent_interval':str(100*(old-current)/old),
 'pinned_gap_closure_percent_interval':str(100*(old-current)/(old-density)),
 'scope':'Formula sourced from published table; arithmetic independently checked. Packing coordinates/construction not replayed; no exhaustive current-best claim.'}
 rows.append(row)
report={'pass':True,'precision_bits':256,'center_density_convention':'Number of unit-radius spheres per unit volume; multiply by unit-ball volume for packing density.','rows':rows}
(D/'KNOWN_LOWER_BOUNDS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
