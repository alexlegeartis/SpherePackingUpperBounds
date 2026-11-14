"""Independent exact period-basis verification of the new57–59D sections."""
from pathlib import Path
import numpy as np,sympy as sp,json,re,hashlib
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;SRC=P.parent/'dual64'
s=(SRC/'CQ64b.html').read_text();body=re.search(r'<a NAME="GRAM".*?</a><br>(.*?)<a NAME="DIVISORS"',s,re.S).group(1)
nums=list(map(int,re.findall(r'-?\d+',re.sub(r'<[^>]*>',' ',body))));assert nums[:2]==[64,0] and len(nums)==2082
G=sp.zeros(64);pos=2
for i in range(64):
 for j in range(i+1):G[i,j]=G[j,i]=nums[pos];pos+=1
assert G.det()==3**32;Gi=G.inv()
def gcdex(a,b):
 oldr,r=a,b;olds,t=1,0;oldt,u=0,1
 while r:
  q=oldr//r;oldr,r=r,oldr-q*r;olds,t=t,olds-q*t;oldt,u=u,oldt-q*u
 if oldr<0:oldr,olds,oldt=-oldr,-olds,-oldt
 return oldr,olds,oldt
def integer_kernel(A):
 B=A.copy();k,n=A.shape;U=sp.eye(n);L=sp.eye(k)
 for i in range(k):
  for j in range(i+1,n):
   a,b=int(B[i,i]),int(B[i,j])
   if b==0:continue
   g,s,t=gcdex(a,b);assert g>0
   ci,cj=B[:,i],B[:,j];ui,uj=U[:,i],U[:,j]
   B[:,i]=s*ci+t*cj;B[:,j]=(-b//g)*ci+(a//g)*cj
   U[:,i]=s*ui+t*uj;U[:,j]=(-b//g)*ui+(a//g)*uj
  if B[i,i]==-1:B[:,i]=-B[:,i];U[:,i]=-U[:,i]
  assert B[i,i]==1
  for j in range(i+1,k):
   a=B[j,i];B[j,:]=B[j,:]-a*B[i,:];L[j,:]=L[j,:]-a*L[i,:]
 assert B==sp.eye(k).row_join(sp.zeros(k,n-k));assert L*A*U==B
 assert abs(U.det())==1 and abs(L.det())==1
 N=U[:,k:].T;assert N*A.T==sp.zeros(n-k,k)
 return N,U,L
out=[]
for name,expected in [('D5',128),('E6',192),('E7',256)]:
 f=P/f'ne64_section_{name}.npz';d=np.load(f);A=sp.Matrix(d['A']);k=A.rows;n=64-k
 assert sp.Matrix(d['G'])==G and sp.Matrix(d['H'])==3*Gi
 K=A*Gi*A.T;assert K==sp.Matrix(d['Knum'])/3 and K.det()==expected
 snf=smith_normal_form(A,domain=sp.ZZ);assert all(abs(snf[i,i])==1 for i in range(k))
 N,U,L=integer_kernel(A);period=N*G*N.T;det=period.det();assert det==3**32*expected
 density_sq=sp.Rational(3**n,det)
 np.savez_compressed(P/f'ne64_section_{name}_period.npz',period=np.array(period,dtype=np.int64),section_basis=np.array(N,dtype=np.int64),column_unimodular=np.array(U,dtype=np.int64),row_unimodular=np.array(L,dtype=np.int64))
 row={'name':name,'dimension':n,'normal_gram_determinant':expected,'section_period_determinant':str(det),'section_kernel_complete':True,'guaranteed_squared_minimum':12,'center_density_squared':str(density_sq),'center_density':float(sp.sqrt(density_sq)),'source_certificate_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'period_file':f'ne64_section_{name}_period.npz'};out.append(row);print(json.dumps(row),flush=True)
report={'all_passed':True,'ambient_minimum_input':12,'minimum_scope':'Inherits published Ne64 minimum12. No new ambient minimum proof.','source_gram_identity_verified':True,'global_record_status':'Not established by this verification.','sections':out}
(P/'ne64_sections_verified.json').write_text(json.dumps(report,indent=2));print('ALL PASS',flush=True)
