#!/usr/bin/env python3
"""A finite, exactly bounded propagation audit for candidate dimensions57--71.

This is NOT an exhaustive record table. Inputs: published Schuett56, Ne64,
Gamma72, low-dimensional lattice/packing densities; exact Gamma72 sections;
and the separately proved Ne64 65D discriminant-count extension.
Families: direct products, torus coloring by K/2K, proper coloring of A_k,
and averaged-ball projections. Vardy's printed rounded numbers stay separate.
All principal computed baseline density SQUARES are rational lower bounds.
"""
from pathlib import Path
import json,math
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
 rows.append({'dimension':n,'name':name,'family':family,'center_density_squared':str(d2),'center_density':float(sp.sqrt(d2)),**kw})
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
# This 65D existential lattice lower bound is independently proved by agents
# higher_parent and code_routes; it is stronger than the older Gamma72 65D row.
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
# Published older rounded values are NOT exact benchmark formulas.
vardy_log={57:15.37,58:16.46,59:17.75,60:19.04,61:20.54,62:22.13,63:23.92}
vardy={n:{'printed_log2_lower_bound':v,'two_to_printed_value':2**v,'scope':'Rounded lower bound printed in Vardy1999, not an exact density identity and not a current world-record claim.'} for n,v in vardy_log.items()}
# Current candidate comparisons use exact certificate values already on disk.
candidates={}
for a in json.loads((P/'ne64_lowrank_sections.json').read_text()):
 candidates[a['dimension']]={'name':'Ne64 bare '+a['name']+' section','d2':sp.sympify(a['center_density_squared'])}
for a in json.loads((R/'dual64'/'constructed.json').read_text()):
 candidates[a['dimension']]={'name':a['label'],'d2':sp.sympify(a['delta_squared'])}
for a in json.loads((R/'higher_parent'/'gamma72_verified.json').read_text())['certificates']:
 n=a['dimension'];d2=sp.sympify(a['center_density_squared_exact'])
 if n not in candidates or d2>candidates[n]['d2']:candidates[n]={'name':a['certificate_file'],'d2':d2}
ne65=json.loads((R/'higher_parent'/'ne65_verified.json').read_text())
assert ne65['status']=='verified_minimum_12_lattice' and ne65['coset_minimum_certified']
assert sp.sympify(ne65['center_density_squared'])==Q(3**34,10)
candidates[65]={'name':'explicit Ne65 (certified coset minimum26/3)','d2':sp.sympify(ne65['center_density_squared'])}
best={};table=[]
for n in range(57,72):
 pool=[a for a in rows if a['dimension']==n]
 b=max(pool,key=lambda a:sp.sympify(a['center_density_squared']));best[n]=b
 c=candidates.get(n);cr=None
 if c:
  cr={'name':c['name'],'center_density_squared':str(c['d2']),'center_density':float(sp.sqrt(c['d2'])),'ratio_to_this_finite_envelope':float(sp.sqrt(c['d2']/sp.sympify(b['center_density_squared']))),'beats_this_finite_envelope':bool(c['d2']>sp.sympify(b['center_density_squared']))}
 table.append({'dimension':n,'best_finite_family_benchmark':b,'candidate':cr,'Vardy1999_rounded':vardy.get(n)})
report={'scope':'Finite propagated lower-bound audit, NOT an exhaustive SOTA/priority determination. It does not bound all sections, laminations, nonlattice packings or constructions from other parents.','published_parent_sources':['https://cohn.mit.edu/sphere-packing/','https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64b.html','https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/Ext72.pdf'],'small_dimension_source':'https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html','Vardy_source':'https://citeseerx.ist.psu.edu/document?doi=c9085ab03f7ab8f5154d3c38fe2be7ab6316fac4&repid=rep1&type=pdf','pi_bounds_certified_by_Machin':True,'minimum_inputs_inherited':True,'all_torus_volume_inequalities_checked_exactly':True,'rows':table,'all_finite_family_bounds':rows}
(P/'propagated_benchmarks.json').write_text(json.dumps(report,indent=2)+'\n')
lines=['# Propagated benchmark audit, dimensions57–71','',report['scope'],'','All densities below are center densities. Principal computed bounds have exact rational squared densities in `propagated_benchmarks.json`. The torus bounds use strict rational volume inequalities, not floating hole searches.','', '| Dimension | Strongest bound in audited finite families | Construction | Current candidate | Vardy1999 printed log₂δ |','|---|---:|---|---:|---:|']
for a in table:
 n=a['dimension'];b=a['best_finite_family_benchmark'];c=a['candidate'];v=a['Vardy1999_rounded']
 lines.append(f"| {n} | {b['center_density']:,.6f} | {b['name']} | {c['center_density']:,.6f} | {v['printed_log2_lower_bound']:.2f}.. |" if c and v else f"| {n} | {b['center_density']:,.6f} | {b['name']} | {c['center_density']:,.6f} | — |" if c else f"| {n} | {b['center_density']:,.6f} | {b['name']} | — | — |")
lines += ['', '## Consequences', '', '- The present Ne64 bare sections in57,58,59 are below explicit propagated bounds. The57D comparison follows from the Schütt56 torus extension; the58D and59D comparisons follow from newly certified primitive Γ72 sections.', '- The63D antipode is below the older printed Vardy1999 lower bound. Its density is therefore not a record candidate.', '- The older Γ72 65D antipode is below the independent Ne64 discriminant-count extension. The candidate column now uses the separately certified explicit Ne65 lattice with center density3^17/√10; this new construction is not counted as a prior baseline.', '- Candidates in60–62 and65–71 exceed these audited finite families. This statement does not establish current SOTA, because these families and the available published tables are not exhaustive.', '', '## Families included', '', '1. Direct products of the published56D and64D packings with listed small-dimensional packings, including the nonlattice10D,11D and13D examples.', '2. Greedy separated offsets on the parent torus, combined with a small-dimensional lattice colored modulo2. The low-dimensional lattice minimum, not a nonlattice density, is used for this family.', '3. The analogous A_k coloring with k+1 colors. Its same-color squared minimum ratio is4 for k=1,3 for k=2,and at least2 otherwise.', '4. Averaged projection through a ball window from64D and72D.', '5. Primitive Γ72 sections along explicitly embedded laminated normals in its √2Leech block.', '6. Direct products of the resulting bounds with listed small-dimensional packings.', '7. The separately proved Ne64 extension with a3-torsion dual coset of minimum at least8.', '', '## Formulas and proof scope', '', 'For a lattice L normalized to minimum4 and center densityδ, greedy torus selection gives s offsets at pairwise distance at leastR whenever (s−1)δ V_n R^n<1. If the vertical min1 latticeK has an s-coloring with same-color squared minimum at leastλ, scale it by t with t²=max(4−R²,4/λ). The packing consisting of L+h_color over each vertical point has center density δ/[t^k covol(K)]. The certificate chooses rational R² and checks the union-volume inequality using π<355/113. Full derivation and an independent57/58 check are in `../code_routes/TORUS_LAYER_BENCHMARKS.md`.', '', 'An averaged k-dimensional ball-window projection of an N-dimensional packing gives δ_(N−k)≥δ_N V_k(k/N)^(k/2)((N−k)/N)^((N−k)/2). These values are evaluated below using the certified lower boundπ>333/106. The formulas are lower bounds, not upper limitations on what projection can achieve.', '', 'For a primitive normal GramK inside unimodular Γ72, the orthogonal section has determinant detK and inherited minimum8. Thus δ²=2^(72−rankK)/detK. The rank13 and14 determinants are2^23 and3·2^22, yielding δ59=262144 and δ58=2^18/√3. Source identity, Smith saturation and explicit complete section-period bases are checked in `verify_gamma72_leech_sections.py`.', '', 'The Vardy values are kept as the printed rounded log₂ lower bounds. No inferred exact Mordell–Weil determinant sequence is used. The catalogue MW64 Gram has determinant2^52, but the attempted identification of its sections with those printed values has not been established.', '', '## Replay', '', '```bash','OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py','OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/propagated_benchmarks.py','```','']
(P/'PROPAGATED_BENCHMARKS.md').write_text('\n'.join(lines))
for a in table:
 b=a['best_finite_family_benchmark'];c=a['candidate'];print(a['dimension'],f"{b['center_density']:.9g}",b['name'],'candidate',None if c is None else round(c['center_density'],6),None if c is None else c['beats_this_finite_envelope'])
