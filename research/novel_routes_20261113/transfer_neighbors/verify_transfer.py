"""Exact finite checks; deliberately excludes floating CVP completeness claims."""
from pathlib import Path
import json,hashlib
from fractions import Fraction
import numpy as np,sympy as sp
P=Path(__file__).resolve().parent
rows=[]
for i in range(7):
 D=np.load(P/f'Q32_{i}_bkz.npz');G=D['G'];H=D['H'];T=D['T'];S=sp.Matrix(G)
 assert np.array_equal(G,np.loadtxt(P/f'Q32_{i}.txt',dtype=np.int64))
 assert np.array_equal(T@G@T.T,H);detT=sp.Matrix(T).det();assert abs(detT)==1
 assert S.det()==65536
 L,d=S.LDLdecomposition(hermitian=False);assert all(d[j,j]>0 for j in range(32))
 rows.append({'index':i,'determinant':65536,'basis_determinant':int(detT),'positive_definite_exact':True})
reject=[]
for r in json.loads((P/'capped_candidate_rejections.json').read_text()):
 G=np.loadtxt(P/f"Q32_{r['index']}.txt",dtype=np.int64);w=np.array(r['rejecting_vector_original_basis'],dtype=np.int64)
 a=np.array([(r['parity']>>i)&1 for i in range(32)])
 assert np.all((w-a)%2==0)
 q=int(w@G@w);assert q==r['exact_rejecting_norm'] and q<=16
 reject.append({'index':r['index'],'exact_rejecting_squared_hole_distance':str(Fraction(q,4))})
delta36=Fraction(16777216,3176523)
bound35sq=delta36**2*Fraction(4,9)*Fraction(35,36)**35
assert bound35sq<8
report={'all_passed':True,'lattices':rows,'false_capped_candidate_rejections':reject,'35d_non_singleton_antipode_upper_bound_squared':str(bound35sq),'35d_non_singleton_antipode_upper_bound':float(bound35sq)**.5,'35d_antipode_scope':'This bounds the standard inherited-minimum mu-beta guarantee, not an independently improved actual minimum of the projected point set.','35d_record_squared':8,'33d_two_periodic_threshold_squared_hole_distance':4,'new_density_record':False,'inherited_inputs':['Catalogue minimum6 for the seven32D parents.','Previously committed36D density.'],'scope':'Exact determinant/congruence/rejecting-witness/analytic-threshold checks. No claim that a floating CVP enumeration is complete.'}
(P/'exact_transfer_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
