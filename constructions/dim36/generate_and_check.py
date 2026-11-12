"""Exact construction, positive definiteness, normalization and verifier inputs."""
from pathlib import Path
import sympy as s
import json
P=Path(__file__).resolve().parent
read=lambda n:s.Matrix([[int(x)for x in r.split()]for r in (P/'data'/n).read_text().splitlines()])
G4,D,U=map(read,['G4.txt','D.txt','U.txt'])
B=5*G4+5*D
assert B==read('Q_original.txt') and B==B.T and abs(U.det())==1
assert all(B[i,i]%2==0 for i in range(36)) and B[0,0]==160
Q=U*B*U.T
L,H=Q.LDLdecomposition(hermitian=False)
assert L*H*L.T==Q and all(H[i,i]>0 for i in range(36))
detphysical=Q.det()/20**36
assert detphysical==s.Rational(10090298369529,4096)
center=2**18/s.sqrt(detphysical)
ratio=s.cancel(center/(s.Rational(2**18,3**10)))
assert center==s.Rational(16777216,3176523) and ratio==s.Rational(139968,117649)
(P/'Q_reduced.txt').write_text('\n'.join(' '.join(str(Q[i,j])for j in range(36))for i in range(36))+'\n')
rows=[]
for i in range(36):
 for x in [H[i,i]]+[L[j,i]for j in range(36)]:
  lo=int(s.floor(x*2**40));hi=int(s.ceiling(x*2**40))
  assert max(abs(lo),abs(hi))<2**53 and s.Rational(lo,2**40)<=x<=s.Rational(hi,2**40)
  rows.append(f'{lo} {hi}')
(P/'ldl.txt').write_text('\n'.join(rows)+'\n')
assert (P/'ldl.txt').read_bytes()==(P/'provenance/certify_ct36_minimum/ldl_quarter.txt').read_bytes()
assert Q==s.Matrix([[int(x)for x in r.split()]for r in (P/'provenance/certify_ct36_minimum/Qquarter.txt').read_text().splitlines()])
print(json.dumps({'dimension':36,'Gram_divisor':20,'minimum_witness':160,'determinant':str(detphysical),'center_density':str(center),'KP36_ratio':str(ratio),'percent_gain':str(s.N(100*(ratio-1),14)),'sympy_version':s.__version__,'exact_checks':'PASS'},indent=2))
