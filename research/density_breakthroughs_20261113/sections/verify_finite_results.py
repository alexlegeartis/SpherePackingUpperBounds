"""Replay exact finite claims; deliberately does not certify floating SVP runs."""
from pathlib import Path
from fractions import Fraction
from math import prod,sqrt
import hashlib,json
import numpy as np
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
out={'scope':'exact finite witnesses, determinants, indices, and arithmetic only; no floating-shell completeness or numerical local-optimality claim'}
cert=json.load(open(P/'mod2_obstruction_certificate.json'))
qfile=ROOT/'constructions/dim36/data/Q_original.txt'
assert hashlib.sha256(qfile.read_bytes()).hexdigest()==cert['integer_gram_sha256']
Q=np.loadtxt(qfile,dtype=np.int64);V=np.array(cert['vectors_original_coordinates'],dtype=np.int64)
assert len(V)==23 and np.all(np.einsum('ij,ij->i',V@Q,V)==160)
assert np.all((V.T@V)%2==0)
out['mod2_rank2_obstruction']={'exact_contact_count':23,'all_contact_squared_norms':8,'outer_product_sum_zero_mod2':True,'conclusion':'No F2-linear map to F2^2 is nonzero on every contact.'}

targets={11:Fraction(3**37,2**42),13:Fraction(3**35,2**38),14:Fraction(196608)}
normal_results=[]
for filename in ['P48p_best_k11.npz','large_best_k12.npz','P48n_extended_best_k13.npz','P48n_extended_best_k14.npz']:
 d=np.load(P/filename);G=d['G'];B=d['B'];k=len(B);K=B@G@B.T
 stored=d['gram'] if 'gram' in d.files else d['K'];assert np.array_equal(K,stored)
 sm=smith_normal_form(s.Matrix(B),domain=s.ZZ);idx=prod(abs(int(sm[i,i]))for i in range(k))
 determinant=Fraction(int(s.Matrix(K).det()),idx*idx);assert determinant.denominator==1
 row={'file':filename,'normal_rank':k,'dimension':48-k,'saturation_index':idx,'determinant':int(determinant),'center_density_using_min6':1.5**((48-k)/2)/sqrt(float(determinant))}
 if k in targets:
  assert determinant>targets[k]
  row['record_beating_determinant_threshold']=str(targets[k]);row['record_beating']=False
 normal_results.append(row)
out['normal_subspaces']=normal_results

# Finite obstruction to all nearest integer roundings of (3/2) Kappa11.
G=np.array(json.load(open(P/'kappa11.json'))['gram'],dtype=np.int64)
assert int(s.Matrix(G).det())==972
raw=np.loadtxt(P/'kappa11_vectors.txt',dtype=np.int64);W=raw[:,1:]
assert np.all(np.einsum('ij,ij->i',W@G,W)==4)
i,j=np.triu_indices(11,1);odd=G[i,j]%2!=0;i=i[odd];j=j[odd];assert len(i)==14
signs=2*((np.arange(2**14)[:,None]>>np.arange(14))&1)-1
slopes=(W[:,i]*W[:,j])@signs.T
admissible=np.all(slopes>=0,axis=0);assert not admissible.any()
out['kappa11_nearest_roundings']={'odd_off_diagonal_entries':14,'roundings_checked':2**14,'valid_necessary_contact_systems':0,'scope':'nearest integer rounding only, using explicit exact norm-four contacts'}

d=np.load(P/'K11_extension.npz');K=d['K'];adj=d['adj'];det=int(d['det'])
assert int(s.Matrix(K).det())==det==121500
assert np.array_equal(K@adj,det*np.eye(11,dtype=np.int64))
cluster_results=[]
for row in json.load(open(P/'K11_clusters.json')):
 S=np.array(row['cluster'],dtype=np.int64);D=S[:,None,:]-S[None,:,:]
 numer=np.einsum('abi,ij,abj->ab',D,adj,D)
 assert len(np.unique(S,axis=0))==row['size'] and numer.max()<=row['beta_numerator']
 beta=Fraction(row['beta_numerator'],row['beta_denominator']);assert 0<beta<=1
 delta=len(S)*float((Fraction(6)-beta)/4)**18.5/sqrt(det)
 assert abs(delta-row['center_density'])<1e-12
 cluster_results.append({'size':len(S),'diameter_squared_bound':str(beta),'center_density':delta})
out['K11_explicit_cluster_checks']=cluster_results
out['K11_cluster_scope']='Checks the displayed clusters, not optimality over the infinite dual lattice. The finite clique searches used a floating-discovered radius-one pool.'

delta36=Fraction(16777216,3176523)
norm35=delta36**2*Fraction(5,4)**35/16
det34=delta36**2*Fraction(5,4)**34/(4*Fraction(3**33,2**50))
out['minimum10_section_thresholds']={'dimension35_dual_norm_squared':str(norm35),'dimension35_threshold_decimal':float(norm35),'dimension34_normal_determinant':str(det34),'dimension34_threshold_decimal':float(det34),'qualification':'Conditional on primitive normal data and an independently established section minimum of at least10; no such section was found.'}
(P/'finite_verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
