from pathlib import Path
from fractions import Fraction as F
import json,math,numpy as np
from flint import fmpq_mat as M,fmpz_mat as Z,fmpq as Q
P=Path(__file__).resolve().parent;S=P/'input'
a=json.loads((S/'candidate78.json').read_text());p=json.loads((S/'period78.json').read_text());r=json.loads((S/'report78.json').read_text());data=np.load(S/'parent/parent80.npz');print(data.files)
G=M(data['G'].tolist());B=M(a['normal_rows']);W=M(p['kernel_rows']);A=M(p['gram']);rep=M(p['fiber_representatives']);shift=M(p['shifts']);K=B*G.inv()*B.transpose()
assert G==G.transpose() and G.det()==41**12
assert K==M([[Q(162,41),Q(-81,41)],[Q(-81,41),Q(162,41)]])
assert W.nrows()==78 and W.ncols()==80 and W*B.transpose()==M(78,2)
assert rep*B.transpose()==M([[0,0],[1,0],[0,1]])
# This explicit unimodular basis proves saturation and exhausts every integer fiber.
U=Z([p['fiber_representatives'][1],p['fiber_representatives'][2]]+p['kernel_rows']);assert abs(U.det())==1
assert W*G*W.transpose()==A and A.det()==G.det()*K.det()==3**9*41**10
labels=M([[0,0],[1,0],[0,1]]);Y=rep-labels*K.inv()*B*G.inv();assert shift*W==Y and Y*B.transpose()==M(3,2)
for i in range(3):
 for j in range(i):
  v=M([[labels[i,k]-labels[j,k]for k in range(2)]])
  assert (v*K.inv()*v.transpose())[0,0]==Q(82,243)
contact=M([r['minimum_attaining_parent_vector']]);assert (contact*G*contact.transpose())[0,0]==14 and contact*B.transpose()==M([[1,0]])
assert list(contact.tolist()[0])==list(data['section_coordinates'][r['minimum_attaining_parent_source_index']])
v=contact-M([[1,0]])*K.inv()*B*G.inv();assert (v*G*v.transpose())[0,0]==Q(3320,243)
delta2=9*F(830,243)**78/F(str(A.det()));assert delta2==F(a['center_density_squared'])
def at(q,N):return sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
pi=16*at(5,201)-4*at(239,200);floor=F('1.41519803035180e-16');assert floor**2<delta2*pi**78/math.factorial(39)**2
out={'status':'PASS exact independent projection/coarea/saturation/shifts/contact/density; parent minimum14 external geometric premise','dimension':78,'floor':str(floor),'center_density_squared':str(delta2),'external_center_density_squared':str(F(2**74,3))}
(P/'result.json').write_text(json.dumps(out,indent=2));print('PASS independent78 geometry and density')
