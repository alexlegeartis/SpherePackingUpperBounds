"""Exact antipode certificates in dimensions 66--71.

Input theorem: the catalogue Gram of Nebe's Gamma72 has minimum 8.
We verify all remaining arithmetic, including source identity, det=1,
positive definiteness, saturation, finite diameters and density ratios.
This is not a world-record/priority certificate.
"""
from pathlib import Path
from fractions import Fraction
import re,json,math,hashlib,time
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent
s=(P/'neb72.html').read_text();t=s[s.index('<a NAME="GRAM"'):s.index('<a NAME="DET"')]
nums=[int(x) for x in re.findall(r'-?\d+',re.sub('<[^>]*>',' ',t))]
assert nums[0]==72 and len(nums)==1+72**2
G=sp.Matrix(72,72,nums[1:]);assert G==G.T
start=time.time();L,D=G.LDLdecomposition(hermitian=False)
assert all(D[i,i]>0 for i in range(72)) and L*D*L.T==G
assert sp.prod(D[i,i] for i in range(72))==1
assert all(G[i,i]%2==0 for i in range(72))
print('Ambient checks passed',time.time()-start,flush=True)
# Exact squared maximal bare-section bound using inherited minimum 8 using known Hermite constants
# in ranks <=8. These apply to every saturated min8 normal lattice.
bare2={1:sp.Integer(2)**68,2:sp.Rational(2**66,3),3:sp.Integer(2)**61,
       4:sp.Integer(2)**58,5:sp.Integer(2)**55,6:sp.Rational(2**54,3),
       7:sp.Integer(2)**50,8:sp.Integer(2)**48}
out=[]
for f in sorted(list(P.glob('gamma72_simplex_rank*.npz'))+list(P.glob('gamma72_root_cluster_rank*.npz'))):
 data=np.load(f);B=sp.Matrix(data['B'].tolist());K=sp.Matrix(data['K'].tolist());S=sp.Matrix(data['S'].tolist());k=B.rows;n=72-k
 assert sp.Matrix(data['G'].tolist())==G and B*G*B.T==K
 assert B.rank()==k
 smith=smith_normal_form(B,domain=sp.ZZ);assert all(abs(smith[i,i])==1 for i in range(k))
 H=K.inv();det=K.det();assert det>0
 beta=max(((S[i,:]-S[j,:])*H*(S[i,:]-S[j,:]).T)[0] for i in range(S.rows) for j in range(S.rows))
 assert beta<8
 # Distinct cosets modulo K are automatic from beta<8 and inherited min8.
 for i in range(S.rows):
  for j in range(i):
   v=(S[i,:]-S[j,:])*H
   assert not all(x.q==1 for x in v)
 d2=sp.Rational(S.rows**2,det)*((8-beta)/4)**n
 ratio2=sp.factor(d2/bare2[k]);assert ratio2>1
 row={'certificate_file':f.name,'dimension':n,'normal_rank':k,'normal_gram':K.tolist(),'normal_determinant':str(det),
      'normal_coefficients':B.tolist(),'saturation_invariants':[1]*k,
      'dual_labels':S.tolist(),'cosets':S.rows,'diameter_squared':str(beta),
      'packing_minimum_squared_at_least':str(8-beta),
      'center_density_squared_exact':str(d2),'center_density_approx':float(sp.sqrt(d2)),
      'inherited_minimum8_section_ceiling_squared_exact':str(bare2[k]),
      'ratio_to_inherited_minimum8_section_ceiling_squared_exact':str(ratio2),
      'ratio_to_inherited_minimum8_section_ceiling_approx':float(sp.sqrt(ratio2)),
      'world_record_established':False,
      'inherited_input':'Published minimum 8 of Nebe Gamma72; not re-enumerated here.'}
 out.append(row);print(n,'det',det,'beta',beta,'density',row['center_density_approx'],'ratio',row['ratio_to_inherited_minimum8_section_ceiling_approx'],flush=True)
report={'ambient_dimension':72,'ambient_determinant':1,'ambient_positive_definite_exact':True,
        'ambient_even_integral_exact':True,'ambient_source_sha256':hashlib.sha256((P/'neb72.html').read_bytes()).hexdigest(),
        'ambient_gram_sha256':hashlib.sha256(np.array(G.tolist(),dtype='<i8').tobytes()).hexdigest(),
        'inherited_minimum':8,'all_checked':True,'world_record_established':False,'comparison_scope':'Bare section density lower bounds computed from parent minimum8, not arbitrary sections with higher actual minimum.', 'certificates':out}
(P/'gamma72_verified.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
print('PASS',len(out),'certificates',time.time()-start,flush=True)
