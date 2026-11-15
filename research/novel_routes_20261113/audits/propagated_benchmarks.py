#!/usr/bin/env python3
"""Exact internal reference arithmetic and explicitly attributed published inputs.

This is NOT an exhaustive record table. Inputs: published Schuett56, Ne64,
Gamma72, low-dimensional lattice/packing densities; exact Gamma72 sections;
and the separately proved Ne64 65D discriminant-count extension.
Families: direct products, torus coloring by K/2K, proper coloring of A_k,
and averaged-ball projections. Vardy's printed rounded numbers stay separate.
All principal computed baseline density SQUARES are rational lower bounds.
"""
from pathlib import Path
import hashlib,json,math
import sympy as sp
import mpmath as mp
P=Path(__file__).resolve().parent;R=P.parent
mp.mp.dps=90
Q=sp.Rational
# Machin alternating-series bounds certify the two rational pi bounds used below.
def atan_bounds(q,n=16):
 terms=[(-1)**j*Q(1,(2*j+1)*q**(2*j+1)) for j in range(n)]
 val=sum(terms);nxt=(-1)**n*Q(1,(2*n+1)*q**(2*n+1))
 return min(val,val+nxt),max(val,val+nxt)
a,b=atan_bounds(5);c,d=atan_bounds(239)
assert Q(333,106)<16*a-4*d and 16*b-4*c<Q(355,113)
pi_lo=Q(333,106);pi_hi=Q(355,113)
seed={56:Q(5**28,2**28*29**4),64:sp.Integer(3)**16,72:sp.Integer(2)**36}
lat2={0:1,1:Q(1,4),2:Q(1,12),3:Q(1,32),4:Q(1,64),5:Q(1,128),6:Q(1,192),7:Q(1,256),8:Q(1,256),9:Q(1,512),10:Q(1,768),11:Q(1,972),12:Q(1,729),13:Q(1,972),14:Q(1,768),15:Q(1,512)}
pack2=dict(lat2);pack2.update({10:Q(25,16384),11:Q(81,65536),13:Q(81,65536)})
rows=[]
def add(n,name,family,d2,**kw):
 d2=sp.cancel(d2)
 assert d2.is_Rational and d2>0
 rows.append({'dimension':n,'name':name,'family':family,'center_density_squared':str(d2),'center_density':float(sp.sqrt(d2)), 'origin': 'published_parent_input' if family=='published_input' else 'computed_in_this_repository', **kw})
for n,d in seed.items():
 if n<=71:add(n,f'published {n}D input','published_input',d*d)
for n,d in seed.items():
 for k in range(1,16):
  if n+k<=71:add(n+k,f'{n}D × best listed {k}D packing','direct_product',d*d*pack2[k],parent=n,added_dimension=k)
# Certified hole radii: evaluate a candidate with high precision, then verify
# the strict union-volume inequality by exact rational arithmetic.
def radius(n,delta,s):
 Vhi=pi_hi**(n//2)/math.factorial(n//2)
 q=1/((s-1)*delta*Vhi)
 raw=mp.power(mp.mpf(int(sp.numer(q)))/int(sp.denom(q)),mp.mpf(2)/n)
 r2=Q(int(mp.floor(raw*10**6)),10**6)
 mass=(s-1)*delta*Vhi*r2**(n//2)
 assert 0<mass<1
 return r2,mass
for n in [56,64]:
 delta=seed[n]
 for k in range(1,min(15,71-n)+1):
  s=2**k;r2,mass=radius(n,delta,s);t2=max(4-r2,Q(1))
  d2=delta**2*4**k*lat2[k]/t2**k
  add(n+k,f'{n}D torus offsets, {k}D lattice modulo 2','torus_binary',d2,parent=n,added_dimension=k,colors=s,horizontal_distance_squared=str(r2),vertical_scale_squared=str(t2),union_volume_upper=str(mass),same_color_minimum_ratio=4)
  # A_k root coloring by sum i*z_i modulo k+1 has no monochromatic root.
  # For k>=3 its same-color min ratio is at least2; k1 and k2 give4 and3.
  s=k+1;lam={1:4,2:3}.get(k,2);r2,mass=radius(n,delta,s);t2=max(4-r2,Q(4,lam))
  a2=Q(1,2**k*(k+1));d2=delta**2*4**k*a2/t2**k
  add(n+k,f'{n}D torus offsets, A{k} with {s} colors','torus_Ak',d2,parent=n,added_dimension=k,colors=s,horizontal_distance_squared=str(r2),vertical_scale_squared=str(t2),union_volume_upper=str(mass),same_color_minimum_ratio=lam)
# Generic projection from a k-ball window: V_k(k/N)^(k/2)(n/N)^(n/2).
# Evaluate V_k with a certified LOWER bound for pi, making d2 rational.
def ball_volume_lower(k):
 if k%2==0:return pi_lo**(k//2)/math.factorial(k//2)
 m=(k-1)//2
 return Q(2**(m+1),math.prod(range(1,2*m+2,2)))*pi_lo**m
for N in [64,72]:
 for n in range(57,min(72,N)):
  k=N-n;v=ball_volume_lower(k)
  d2=seed[N]**2*v*v*Q(k,N)**k*Q(n,N)**n
  add(n,f'{N}D averaged ball-window projection','generic_projection',d2,parent=N,codimension=k)
sections=json.loads((P/'gamma72_leech_sections.json').read_text())
for sn,row in sections.items():
 n=int(sn)
 if 57<=n<=71:add(n,f'Gamma72 laminated normal rank {72-n}','Gamma72_laminated_section',sp.sympify(row['center_density_sq']),certificate=f'gamma72_leech_section_rank{72-n}.npz')
# Classical kappa13 is denser than the laminated normal used above. Its
# published catalogue embedding is checked independently of our search.
literature59_path=P.parents[2]/'constructions/dim59/literature_section_comparison.json'
literature59=json.loads(literature59_path.read_text())
assert literature59['status']=='verified_classical_section_derived_from_published_inputs'
assert literature59['kappa13_embedding_primitive'] and not literature59['is_reported_59D_result']
assert sp.Rational(literature59['center_density_squared_exact'])==Q(2**44,3**5)
add(59,'Gamma72 perpendicular to sqrt(2)*kappa13 (derived from N12 + K13)',
 'Gamma72_classical_kappa13_section',Q(2**44,3**5),
 certificate='../../../constructions/dim59/literature_section_comparison.json',
 reference_kind='derived_from_published_inputs_not_a_quoted_59D_record')
# This existential extension was proved in this repository, in higher_parent
# and code_routes. It is not a published same-dimensional comparison result.
add(65,'Ne64 discriminant-count extension (coset minimum≥8)','Ne64_existential_extension',Q(3**33,4),proof='../code_routes/NE64_THETA_COUNT.md')
# Close under one product with a known packing of dimension<=15. This also
# includes products of torus extensions and primitive Gamma72 sections.
initial=list(rows)
for a in initial:
 if a['family']=='direct_product':continue
 for k in range(1,16):
  n=a['dimension']+k
  if n>71:break
  add(n,a['name']+f' × {k}D packing','propagated_product',sp.sympify(a['center_density_squared'])*pack2[k],parent_construction=a['name'],added_dimension=k)
# Provenance is explicit: the following are OUR comparison constructions.
# Only Vardy's historical table below is a published same-dimensional comparison.
from benchmark_report import render
ROOT = P.parents[2]
sources = json.loads((P/'benchmark_sources.json').read_text())
assert {'V99','N98','N99','N12','S19','CS96','TABLE'} <= sources.keys()
for source in sources.values():
 if 'archived_source' in source: assert (P/source['archived_source']).is_file()

vardy_log = {57:'15.37…',58:'16.46…',59:'17.75…',60:'19.04…',61:'20.54…',62:'22.13…',63:'23.92…'}
published = {n:{'source_id':'V99', 'construction':'Elkies section E_'+str(n),
 'locator':sources['V99']['locator'], 'printed_log2_center_density':value,
 'exact_center_density_squared':None,
 'scope':'Historical printed lower bound, not an exact density identity or a current record.'}
 for n,value in vardy_log.items()}
published[64] = {'source_ids':['N98','N99'], 'construction':'Ne64',
 'exact_center_density_squared':str(3**32), 'scope':'Published 64D parent, also used as a failed-control comparator.'}

candidates = {}
def candidate(n,name,d2,path):
 assert (P/path).is_file(), path
 candidates[n] = {'name':name,'d2':sp.Rational(d2),'certificate_path':path}
for a in json.loads((P/'ne64_lowrank_sections.json').read_text()):
 candidate(a['dimension'],'Ne64 bare '+a['name']+' section',a['center_density_squared'],f"ne64_section_{a['name']}.npz")
for a in json.loads((R/'dual64'/'constructed.json').read_text()):
 candidate(a['dimension'],f"Ne64 antipode, {a['cosets']} translates",a['delta_squared'],f"../dual64/{a['label']}.npz")
for a in json.loads((R/'higher_parent'/'gamma72_verified.json').read_text())['certificates']:
 n=a['dimension']; d2=sp.Rational(a['center_density_squared_exact'])
 if n not in candidates or d2>candidates[n]['d2']:
  candidate(n,f"Gamma72 antipode, {a['cosets']} translates",d2,'../higher_parent/'+a['certificate_file'])
ne65 = json.loads((R/'higher_parent'/'ne65_verified.json').read_text())
assert ne65['status']=='verified_minimum_12_lattice' and ne65['coset_minimum_certified']
assert sp.Rational(ne65['center_density_squared'])==Q(3**34,10)
candidate(65,'Explicit Ne65, certified depth 26/3',ne65['center_density_squared'],'../higher_parent/ne65_verified.json')
# Include the result added on November 15: never regenerate a stale 59D row.
ne59path = ROOT/'constructions/dim59/result.json'
ne59 = json.loads(ne59path.read_text())
ne59audit = json.loads((ne59path.parent/'independent_audit.json').read_text())
assert ne59['status']=='verified_periodic_packing' and ne59['all_checks_exact']
assert ne59audit['all_checks_passed'] and ne59['data_sha256']==ne59audit['data_sha256']
assert ne59['center_density_squared_exact']==ne59audit['center_density_squared_exact']
assert sp.Rational(ne59['center_density_squared_exact'])==Q(3**33*401**59,2**296*5**2*13**61)
candidate(59,'Six-translate Ne64 antipode',ne59['center_density_squared_exact'],'../../../constructions/dim59/result.json')

# Each computed reference names the actual local evidence, not a parent paper.
for a in rows:
 family=a['family']
 if family in ('Gamma72_laminated_section','Gamma72_classical_kappa13_section'):a['evidence_path']=a['certificate']
 elif family=='Ne64_existential_extension':a['evidence_path']=a['proof']
 elif family=='published_input':a['evidence_path']='benchmark_sources.json'
 else:a['evidence_path']='propagated_benchmarks.py'
 assert (P/a['evidence_path']).is_file()
 if family=='Gamma72_classical_kappa13_section':a['parent_source_ids']=['N12','K13']
 elif family=='Gamma72_laminated_section':a['parent_source_ids']=['N12']
 elif a.get('parent')==56:a['parent_source_ids']=['S19']
 elif a.get('parent')==64 or family=='Ne64_existential_extension':a['parent_source_ids']=['N98','N99']
 elif a.get('parent')==72:a['parent_source_ids']=['N12']

best = {}; table = []
for n in range(57,72):
 pool=[a for a in rows if a['dimension']==n]
 b=max(pool,key=lambda a:sp.Rational(a['center_density_squared']));best[n]=b
 c=candidates[n]; ratio2=sp.cancel(c['d2']/sp.Rational(b['center_density_squared']))
 cr={'name':c['name'],'certificate_path':c['certificate_path'],
     'center_density_squared':str(c['d2']),'center_density':float(sp.sqrt(c['d2']))}
 internal={'ratio_squared_exact':str(ratio2),'ratio':float(sp.sqrt(ratio2)),
     'increase_percent':float(100*(sp.sqrt(ratio2)-1)),
     'exceeds_reference':bool(ratio2>1),
     'comparison_kind':'published_parent_control' if b['origin']=='published_parent_input' else 'internal_repository_comparison',
     'is_gain_over_identified_published_record':False}
 if b['family']=='Gamma72_classical_kappa13_section':
  internal['comparison_kind']='derived_from_published_inputs_not_a_quoted_59D_record'
 table.append({'dimension':n,'computed_reference':b,'repository_construction':cr,
     'internal_comparison':internal,'published_reference':published.get(n),
     'current_published_record_identified':False,
     'published_record_improvement_established':False})
assert sp.Rational(best[59]['center_density_squared'])==Q(2**44,3**5)
assert table[2]['repository_construction']['center_density_squared']==ne59['center_density_squared_exact']
external59 = json.loads((ne59path.parent/'external_comparison.json').read_text())
assert external59['status']=='verified_gain_interval_over_reported_guarantee'
assert external59['construction_result_sha256']==hashlib.sha256(ne59path.read_bytes()).hexdigest()
assert not external59['improvement_over_actual_E59_established']
table[2]['published_guarantee_comparison'] = {
 'kind':'verified_interval_over_reported_lower_bound_not_actual_lattice_density',
 'gain_percent_interval':external59['gain_percent_interval'],
 'evidence_path':'../../../constructions/dim59/external_comparison.json',
 'verifier_path':'../../../constructions/dim59/verify_external_comparison.py'}
report={'schema_version':2,
 'scope':'Exact comparisons with our computed reference constructions; published parent inputs and historical same-dimensional values attributed separately. No improvement over the current published record is established by this audit.',
 'sources':sources,'pi_bounds_certified_by_Machin':True,'minimum_inputs_inherited':True,
 'all_torus_volume_inequalities_checked_exactly':True,
 'rows':table,'all_finite_family_bounds':rows}
(P/'propagated_benchmarks.json').write_text(json.dumps(report,indent=2)+'\n')
(P/'PROPAGATED_BENCHMARKS.md').write_text(render(report))
for a in table:
 print(a['dimension'],a['repository_construction']['name'],
       f"{a['repository_construction']['center_density']:.9f}",
       a['internal_comparison']['comparison_kind'],
       'published record improvement:',a['published_record_improvement_established'])
