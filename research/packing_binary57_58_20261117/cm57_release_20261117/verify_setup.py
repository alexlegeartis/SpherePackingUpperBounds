"""Exact setup verification; does not substitute for complete enumeration."""
from pathlib import Path
from fractions import Fraction as F
import json, math, struct
import numpy as np
from flint import fmpz_mat, fmpq_mat
p=Path(__file__).parent
read=lambda f:json.load(open(p/f))
base=fmpz_mat(np.load(p/'parent/parent56.npz')['G'].tolist())
a=read('cvp_basis_actions.json'); T=fmpz_mat(a['T']); Q=fmpz_mat(a['G']); C=fmpz_mat(a['old_to_new'])
assert abs(T.det())==1 and T*base*T.transpose()==Q and base.det()==29**8
old=read('parent/primal_cm.json'); OT=fmpz_mat(old['T']); OQ=fmpz_mat(old['G']); CM=fmpz_mat(old['CM']); M=fmpz_mat(old['M']); I=fmpz_mat([[int(i==j) for j in range(56)] for i in range(56)])
assert C*T==OT and OT*base*OT.transpose()==OQ
assert CM*OQ*CM.transpose()==OQ and M*(I-CM)==3*I and M*OQ*M.transpose()==3*OQ
cl=read('target_class.json'); row=fmpz_mat([cl['representative']]); pre=fmpz_mat([cl['preimage']]); assert all(int(x)%3==0 for x in (pre*M-row).entries())
e=read('random56_37_exact.json'); r=e['target_numerator']; n=len(r)
assert e['G']==a['G'] and n==56 and e['threshold_numerator']==54 and e['denominator']==3
assert [int(x)%3 for x in (row*C).entries()]==r
D=list(map(F,e['diagonal_rational'])); mu=[[F(x) for x in z] for z in e['mu_rational']]
assert all(x>0 for x in D)
for i in range(n):
 for j in range(n):assert sum(mu[i][k]*D[k]*mu[j][k] for k in range(n))==int(Q[i,j])
lines=(p/'random56_37_bits.txt').read_text().splitlines(); assert int(lines[0])==n and int(lines[1])==54 and list(map(int,lines[2].split()))==r
assert [list(map(int,z.split())) for z in lines[3:3+n]]==a['G']
vals=D+[v for z in mu for v in z]
assert len(lines)==3+n+len(vals)
for line,v in zip(lines[3+n:],vals):
 lo,hi=[struct.unpack('>d',struct.pack('>Q',int(x)))[0] for x in line.split()]
 assert math.isfinite(lo) and math.isfinite(hi) and F.from_float(lo)<=v<=F.from_float(hi)
assert all(Q[i,i]%2==0 for i in range(n)) and (fmpz_mat([r])*Q*fmpz_mat([r]).transpose())[0,0]%6==0
contact=read('random37_contact.json'); v=fmpz_mat([contact['residual_numerator']]); w=fmpz_mat([contact['parent_coordinate_numerator']])
assert v*T==w and v*Q*v.transpose()==w*base*w.transpose() and int((v*Q*v.transpose())[0,0])==60
assert all((int(v[0,i])-r[i])%3==0 for i in range(n))
# Machin identity, alternating rational series: odd partial sum is a lower bound.
def atan_lo(q,N=121):return sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N+1)),F(0))
def atan_hi(q,N=120):return sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N+1)),F(0))
pi_lo=16*atan_lo(5)-4*atan_hi(239)
results=[]
for z in read('extension_candidates.json'):
 dim=z['dimension']; A=fmpq_mat(z['gram']); det=F(29**8)*(F(10,3) if dim==57 else F(25,3)); assert F(str(A.det()))==det==F(z['gram_determinant'])
 for i in range(56):
  for j in range(56): assert A[i,j]==base[i,j]
  assert F(str(A[i,56]))==F(int((w*base)[0,i]),3)
  if dim==58:assert A[i,57]==-A[i,56]
 assert A[56,56]==10
 if dim==58:assert A[57,57]==10 and A[56,57]==A[57,56]==-5
 density2=F(5,2)**dim/det; assert density2==F(z['center_density_squared'])
 if dim%2==0: volume2=pi_lo**dim/F(math.factorial(dim//2)**2)
 else:
  k=(dim-1)//2; volume2=F(2**(k+1),math.prod(range(1,dim+1,2)))**2*pi_lo**(2*k)
 floor=F('1.53307111021240e-11' if dim==57 else '5.02419365648403e-12'); assert floor**2<volume2*density2
 results.append({'dimension':dim,'rigorous_density_floor':str(floor),'determinant':str(det)})
print(json.dumps({'status':'EXACT_SETUP_AND_DENSITY_PASS','enumeration_status':'Requires separate complete run and independent audit','results':results}))
