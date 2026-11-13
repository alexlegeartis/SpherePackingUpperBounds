"""Independent exact audit of the four Ne64 antipode constructions.

Does not import the construction/search code. Reads original archived HTML and
finite witnesses, builds actual section-period bases and all coset shifts, and
checks their determinant. Published ambient minimum12 remains an explicit input.
"""
from pathlib import Path
import re,html,json,itertools,hashlib,time,math
import numpy as np,sympy as sp
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;SRC=P.parent/'dual64'
source=(SRC/'CQ64b.html').read_text()
body=re.search(r'<a NAME="GRAM".*?</a><br>(.*?)<a NAME="DIVISORS"',source,re.S).group(1)
nums=list(map(int,re.findall(r'-?\d+',html.unescape(re.sub(r'<[^>]*>',' ',body)))))
assert nums[:2]==[64,0] and len(nums)==2082
G=sp.zeros(64);k=2
for i in range(64):
 for j in range(i+1):G[i,j]=G[j,i]=nums[k];k+=1
assert G.det()==3**32
L,D=G.LDLdecomposition(hermitian=False);assert all(D[i,i]>0 for i in range(64))
Gi=G.inv();H=3*Gi;assert all(x.q==1 for x in H)
records=json.loads((SRC/'constructed.json').read_text());checks=[]
for row in records:
 label=row['label'];F=SRC/f'{label}.npz';d=np.load(F)
 assert sp.Matrix(d['G'])==G and sp.Matrix(d['H'])==H
 A=sp.Matrix(d['A']);S=sp.Matrix(d['S']);Knum=sp.Matrix(d['Knum']);r=A.rows;n=64-r
 assert A*H*A.T==Knum
 K=A*Gi*A.T;assert K==Knum/3;Ki=K.inv()
 # An explicit unimodular minor is a stronger, constructive saturation check.
 unit=None
 for cols in itertools.combinations(range(64),r):
  B=A[:,list(cols)]
  if abs(B.det())==1:unit=list(cols);break
 assert unit is not None
 B=A[:,unit];Bi=B.inv();free=[j for j in range(64)if j not in unit]
 right=sp.zeros(64,r)
 for j,c in enumerate(unit):right[c,:]=Bi[j,:]
 assert A*right==sp.eye(r) and all(x.q==1 for x in right)
 snf=smith_normal_form(A,domain=sp.ZZ);assert all(abs(snf[i,i])==1 for i in range(r))
 # Every row of N is an integer coordinate row in the original primal basis.
 # Since its free-coordinate block is I, these rows generate the whole kernel.
 N=sp.zeros(n,64);left=-A[:,free].T*Bi.T
 for j,c in enumerate(unit):N[:,c]=left[:,j]
 for j,c in enumerate(free):N[j,c]=1
 assert N*A.T==sp.zeros(n,r) and all(x.q==1 for x in N)
 period=N*G*N.T
 period_det=period.det()
 assert period_det==3**32*K.det()
 lifts=S*right.T;assert lifts*A.T==S and all(x.q==1 for x in lifts)
 Y=lifts-S*Ki*A*Gi;assert Y*A.T==sp.zeros(S.rows,r)
 shifts=Y[:,free];assert shifts*N==Y
 beta=sp.Rational(0)
 distances=[]
 for i,j in itertools.combinations(range(S.rows),2):
  ds=S[i,:]-S[j,:];q=(ds*Ki*ds.T)[0];assert q>0
  distances.append(str(q));beta=max(beta,q)
 assert beta<12
 # If two section translates were equal, a nonzero vector of the ambient
 # lattice would lie in the normal space with this squared norm <= beta<12.
 # That contradicts the inherited ambient minimum, so the cosets are distinct.
 mu=12-beta
 delta_sq=sp.Rational(S.rows**2,1)*(mu/4)**n/period_det
 assert n==row['dimension'] and S.rows==row['cosets']
 assert beta==sp.Rational(row['beta']) and K.det()==sp.Rational(row['normal_det'])
 assert delta_sq==sp.Rational(row['delta_squared'])
 assert abs(float(sp.sqrt(delta_sq))-row['delta'])<1e-8
 den=sp.ilcm(*[x.q for x in shifts]);shnum=shifts*den
 np.savez_compressed(P/f'{label}_explicit_period.npz',period=np.array(period,dtype=np.int64),section_basis=np.array(N,dtype=np.int64),ambient_lifts=np.array(lifts,dtype=np.int64),shift_numerators=np.array(shnum,dtype=np.int64),shift_denominator=np.array(int(den),dtype=np.int64),unit_minor_columns=np.array(unit,dtype=np.int64))
 checks.append({'label':label,'dimension':n,'normal_rank':r,'source_gram_identity':True,'dual_identity_H_equals_3Ginv':True,'unit_minor_columns':unit,'unit_minor_determinant':int(B.det()),'integer_label_lifts_verified':True,'section_kernel_basis_verified':True,'normal_gram_determinant':str(K.det()),'section_period_determinant':str(period_det),'period_covolume_formula':'3^16 * sqrt(detK)','cosets':S.rows,'cluster_squared_distances':distances,'beta':str(beta),'coset_distinctness_reason':'All nonzero normal differences have squared norm at most beta<12, so none can be a nonzero ambient lattice vector.','guaranteed_squared_minimum':str(mu),'center_density_squared':str(delta_sq),'center_density':float(sp.sqrt(delta_sq)),'explicit_period_file':f'{label}_explicit_period.npz','shift_denominator':int(den),'certificate_sha256':hashlib.sha256(F.read_bytes()).hexdigest()})
 print(label,'PASS','density',float(sp.sqrt(delta_sq)),flush=True)
out={'all_passed':True,'source':'https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html','source_sha256':hashlib.sha256((SRC/'CQ64b.html').read_bytes()).hexdigest(),'ambient_determinant':str(G.det()),'ambient_minimum_input':12,'minimum_scope':'Published Ne64 minimum12 is inherited, not independently re-proved here.','global_record_status':'Not established by this mathematical audit.','checks':checks}
(P/'ne64_antipodes_audit.json').write_text(json.dumps(out,indent=2));print('ALL PASS',flush=True)
