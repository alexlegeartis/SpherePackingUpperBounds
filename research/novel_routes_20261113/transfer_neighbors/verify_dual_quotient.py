"""Exact verification of the finite L*/L covering witness."""
from pathlib import Path
import argparse,json
from fractions import Fraction
import numpy as np,sympy as sp
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('index',type=int,nargs='?');args=ap.parse_args()
indices=range(7)if args.index is None else [args.index];out=[]
for i in indices:
 f=P/f'dual_quotient_{i}.npz'
 if not f.exists():continue
 s=np.load(f);G=s['G'].astype(np.int64);D=s['D'].astype(np.int64);H=s['H'].astype(np.int64);T=s['T'].astype(np.int64);C=s['C'].astype(np.int64);piv=s['piv'];W=s['witnesses'].astype(np.int64);covered=s['covered']
 assert np.array_equal(G,np.loadtxt(P/f'Q32_{i}.txt',dtype=np.int64));assert sp.Matrix(G).det()==65536
 assert np.array_equal(G@D,2*np.eye(32,dtype=np.int64));assert np.array_equal(T@D@T.T,H);assert abs(sp.Matrix(T).det())==1;assert np.array_equal(C,T@D)
 assert np.all(np.diag(G)%2==0) and np.all(np.diag(D)%2==0)
 labs=np.zeros(len(W),dtype=np.uint16);tags=np.array([sum(int(C[j,k]%2)<<a for a,k in enumerate(piv))for j in range(32)],dtype=np.uint16)
 for j in range(32):labs^=np.where((W[:,j]&1)!=0,tags[j],np.uint16(0))
 assert np.array_equal(labs[covered],np.flatnonzero(covered));norm=np.einsum('bi,ij,bj->b',W,H,W);assert np.all(norm[covered]<=8)
 assert covered.all() and set(np.unique(norm))=={0,6,8}
 E=W[1<<np.arange(16)];B=(E@H@E.T)%2;A=B.copy();rank=0
 for k in range(16):
  ix=np.flatnonzero(A[rank:,k])
  if not len(ix):continue
  j=rank+ix[0];A[[rank,j]]=A[[j,rank]]
  for l in range(16):
   if l!=rank and A[l,k]:A[l]^=A[rank]
  rank+=1
 assert rank==16 and np.all(np.diag(B)==0)
 bits=((np.arange(65536)[:,None]>>np.arange(16))&1).astype(np.uint8)
 qb=(np.einsum('bi,ij,bj->b',E,H,E)//2)%2;qs=np.zeros(65536,dtype=np.uint8)
 for j in range(16):
  qs^=bits[:,j]*np.uint8(qb[j])
  for k in range(j):
   if B[j,k]:qs^=bits[:,j]*bits[:,k]
 assert np.array_equal(qs,(norm//2)%2)
 assert np.all(norm[qs==1]==6) and np.all(norm[qs==0]<=8)
 vals,cts=np.unique(norm[covered],return_counts=True)
 row={'index':i,'classes':65536,'exactly_covered_classes':int(covered.sum()),'all_classes_covered':bool(covered.all()),'witness_squared_norm_distribution':{str(sp.Rational(int(v),2)):int(c)for v,c in zip(vals,cts)},'quadratic_polar_rank_mod2':rank,'q1_cosets_have_norm3_witness':True,'q0_cosets_have_norm_at_most4_witness':True,'max_witness_abs_coordinate':int(np.abs(W).max()),'scope':'Upper bound4 on the minimum norm in each covered dual coset, not a global covering-radius claim.'};out.append(row);print(json.dumps(row),flush=True)
(P/'dual_quotient_verified.json').write_text(json.dumps(out,indent=2))
mu1=Fraction(3*65536,65535);mu2=Fraction(4*65536-3,65535);record_sq=Fraction(3**33,2**50)
endpoints=[(mu1,mu1**32/Fraction(2**50)),(mu2,mu2**33/Fraction(2**66)/(mu2-3)),(Fraction(6),record_sq)]
assert all(value<=record_sq for _,value in endpoints)
theorem={'scope':'Any packing in33D that is a union of full translates of a fixed catalogue32D parent, with translation offsets in its dual lattice and arbitrary real layer heights. The overall minimum may be any positive value, and the layers need not be equally spaced or periodic.','record_squared_center_density':str(record_sq),'parents_verified':[x['index']for x in out],'maximum_pairwise_q_zero_colors':256,'intensity_bounds':['65536/sqrt(mu)','256/sqrt(mu-3) for mu>3','1/sqrt(mu-4) for mu>4'],'breakpoints':[str(mu1),str(mu2)],'endpoint_squared_density_ratios':[{'mu':str(mu),'ratio_float':float(value/record_sq),'ratio_exact':str(value/record_sq)}for mu,value in endpoints],'derivative_sign_for_mu_power33_div_mu_minus_a':'32*mu-33*a; every interior critical point is a minimum','strict_record_improvement_impossible_in_stated_family':True,'new_density_record':False}
(P/'layer33_theorem_checked.json').write_text(json.dumps(theorem,indent=2));print(json.dumps({k:v for k,v in theorem.items()if k!='endpoint_squared_density_ratios'}),flush=True)
