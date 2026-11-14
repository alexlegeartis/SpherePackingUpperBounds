#!/usr/bin/env python3
"""Exact verification of primitive Gamma72 sections from a laminated Leech chain.

The sole minimum input is Nebe's published min(Gamma72)=8 theorem.
The discovery enumeration is neither read nor trusted. Source Gram identity,
positive definiteness, determinant, all normal Grams and saturation are checked.
For ranks13 and14, explicit full integer period bases are also constructed.
"""
from pathlib import Path
import re,json,hashlib
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent
source=P.parent/'higher_parent'/'neb72.html'
s=source.read_text();body=s[s.index('<a NAME="GRAM"'):s.index('<a NAME="DET"')]
nums=list(map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',body))))
assert nums[0]==72 and len(nums)==1+72**2
G=sp.Matrix(72,72,nums[1:]);assert G==G.T
L,D=G.LDLdecomposition(hermitian=False)
assert all(D[i,i]>0 for i in range(72)) and L*D*L.T==G
assert sp.prod(D[i,i] for i in range(72))==1
F=G[:24,:24]/2
assert all(x.q==1 for x in F) and all(F[i,i]%2==0 for i in range(24)) and F.det()==1

def gcdex(a,b):
 r0,r1,s0,s1,t0,t1=a,b,1,0,0,1
 while r1:
  q=r0//r1;r0,r1=r1,r0-q*r1;s0,s1=s1,s0-q*s1;t0,t1=t1,t0-q*t1
 if r0<0:r0,s0,t0=-r0,-s0,-t0
 return r0,s0,t0

def integer_kernel(A):
 B=A.copy();k,n=A.shape;U=sp.eye(n);E=sp.eye(k)
 for i in range(k):
  for j in range(i+1,n):
   a,b=int(B[i,i]),int(B[i,j])
   if not b:continue
   g,s,t=gcdex(a,b)
   ci,cj=B[:,i],B[:,j];ui,uj=U[:,i],U[:,j]
   B[:,i]=s*ci+t*cj;B[:,j]=(-b//g)*ci+(a//g)*cj
   U[:,i]=s*ui+t*uj;U[:,j]=(-b//g)*ui+(a//g)*uj
  if B[i,i]==-1:B[:,i]=-B[:,i];U[:,i]=-U[:,i]
  assert B[i,i]==1
  for j in range(i+1,k):
   a=B[j,i];B[j,:]-=a*B[i,:];E[j,:]-=a*E[i,:]
 assert B==sp.eye(k).row_join(sp.zeros(k,n-k)) and E*A*U==B
 assert abs(U.det())==abs(E.det())==1
 N=U[:,k:].T;assert N*A.T==sp.zeros(n-k,k)
 return N,U,E

out=[]
for k in range(1,17):
 f=P/f'gamma72_leech_section_rank{k}.npz';data=np.load(f)
 B=sp.Matrix(data['B']);K=sp.Matrix(data['K']);n=72-k
 assert sp.Matrix(data['G'])==G and B.shape==(k,72) and B*G*B.T==K
 assert B[:,24:]==sp.zeros(k,48)
 sf=smith_normal_form(B,domain=sp.ZZ);assert all(abs(sf[i,i])==1 for i in range(k))
 det=K.det();d2=sp.Rational(2**n,det)
 row={'file':f.name,'dimension':n,'normal_rank':k,'normal_determinant':str(det),'primitive':True,'minimum_at_least':8,'center_density_squared':str(d2),'center_density':float(sp.sqrt(d2)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
 if k in [13,14]:
  Nu,U,E=integer_kernel(B);N=Nu*G.inv();assert N*G*B.T==sp.zeros(n,k);period=N*G*N.T
  assert period.det()==det
  np.savez_compressed(P/f'gamma72_leech_section_rank{k}_period.npz',section_basis=np.array(N,dtype=str),period=np.array(period,dtype=str),column_unimodular=np.array(U,dtype=str),row_unimodular=np.array(E,dtype=str))
  row['complete_period_basis_checked']=True
 out.append(row);print(json.dumps(row),flush=True)
report={'all_passed':True,'source':str(source.relative_to(P.parent)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'minimum_input':'Published minimum8 theorem for Gamma72; not re-proved here.','scope':'Explicit valid lower bounds. No priority or global optimality claim.','source_gram_identity':True,'ambient_det':1,'leech_block_det':1,'leech_block_minimum_at_least':4,'sections':out}
(P/'gamma72_leech_sections_verified.json').write_text(json.dumps(report,indent=2)+'\n')
print('ALL PASS',flush=True)
