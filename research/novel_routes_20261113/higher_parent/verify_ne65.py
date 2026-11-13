#!/usr/bin/env python3
"""Verify explicit Ne64 extension; require the independent CVP certificate.

Without ../ne64_hole_certificate/summary.json this reports only conditional
algebraic checks. It never converts exploratory node caps into a proof.
"""
from pathlib import Path
import argparse,hashlib,json,re
import numpy as np
import sympy as s
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--certificate',type=Path,default=P.parent/'ne64_hole_certificate'/'summary.json')
args=ap.parse_args()
data=np.load(P/'ne65_candidate.npz')
G=s.Matrix(data['G64'].tolist());A=s.Matrix(data['G65'].tolist())
w=s.Matrix([data['coset_numerator_original'].tolist()]);z=s.Matrix([data['coset_numerator_reduced'].tolist()])
T=s.Matrix(data['T'].tolist());R=s.Matrix(data['R'].tolist())
a=data['coset_original_residue'].tolist();c=data['coset_reduced_residue'].tolist()
assert G.shape==(64,64) and A.shape==(65,65)
# Match the actual archived published parent, not just another stored Gram.
source=(P/'CQ64b.html').read_text()
raw=source[source.index('<a NAME="GRAM"'):source.index('<a NAME="DIVISORS"')]
vals=[int(x) for x in re.findall(r'-?\d+',re.sub('<[^>]*>',' ',raw))]
assert vals[:2]==[64,0] and len(vals)==2082
Gs=s.zeros(64);k=2
for i in range(64):
    for j in range(i+1):Gs[i,j]=Gs[j,i]=vals[k];k+=1
assert G==Gs and G==G.T
L,D=G.LDLdecomposition(hermitian=False)
assert all(D[i,i]>0 for i in range(64))
assert s.prod(D[i,i] for i in range(64))==3**32
assert all(G[i,i]%2==0 for i in range(64))
Ti=T.inv();assert all(x.q==1 for x in Ti)
assert R==T*G*T.T and w==z*T
assert list(w.applyfunc(lambda x:x%3))==a
assert list(z.applyfunc(lambda x:x%3))==c
assert any(a) and any(c)
assert (w*G*w.T)[0]==(z*R*z.T)[0]==78
assert all(x%3==0 for x in w*G)
H=3*G.inv();assert all(x.q==1 for x in H)
assert (w*G*w.T)[0]%18==6
h=w/3
assert (h*G*h.T)[0]==s.Rational(26,3)
assert all(x.q==1 for x in h*G)
assert A[:64,:64]==G
assert A[:64,64]==G*h.T and A[64,:64]==h*G
assert A[64,64]==12 and A==A.T
assert all(x.q==1 for x in A) and all(A[i,i]%2==0 for i in range(65))
schur=A[64,64]-(h*G*h.T)[0]
assert schur==s.Rational(10,3)>0
assert A.det(method='domain-ge')==10*3**31
# In nonzero classes, the certified horizontal bound combines with m^2*t^2.
assert s.Rational(26,3)+s.Rational(10,3)==12
assert 3**2*s.Rational(10,3)==30>12
# Layer0 has minimum12 by the published parent theorem; the new generator
# itself has norm12. Thus the lower and upper minimum claims coincide.
certified=False;certificate_info=None
if args.certificate.exists():
    cert=json.loads(args.certificate.read_text())
    assert cert.get('certified_empty') is True, 'CVP summary must explicitly certify emptiness'
    cp=args.certificate.parent
    setup=json.loads((cp/'setup.json').read_text())
    assert setup['inclusive_integer_norm_bound']==60
    assert setup['congruence_modulus']==3 and setup['class_id']==4
    assert setup['original_residue']==a and setup['reduced_residue']==c
    assert setup['residue_norm_mod_18']==6
    qr=s.Matrix([[int(x) for x in line.split()] for line in (cp/'gram.txt').read_text().splitlines()])
    cr=[int(x) for x in (cp/'coset.txt').read_text().split()]
    assert qr==R and cr==c
    for name,digest in setup['sha256'].items():
        assert hashlib.sha256((cp/name).read_bytes()).hexdigest()==digest
    certified=True
    certificate_info={'path':str(args.certificate),'sha256':hashlib.sha256(args.certificate.read_bytes()).hexdigest()}
report={
    'status':'verified_minimum_12_lattice' if certified else 'conditional_algebra_verified_CVP_certificate_missing',
    'dimension':65,'parent':'Ne64','parent_gram_matches_archived_catalogue':True,
    'parent_minimum_assumption':'Published Ne64 minimum12 theorem; not re-enumerated here',
    'dual_coset_numerator_norm':78,'dual_coset_norm':'26/3',
    'coset_minimum_certified':certified,'coset_minimum_if_certified':'26/3',
    'height_squared':'10/3','gram_even_integral':True,'gram_positive_definite':True,
    'determinant':str(10*3**31),'determinant_formula':'10*3^31',
    'minimum':12 if certified else None,
    'layer_bounds':{'zero_layer':'12 from published parent minimum','nonzero_mod3':'26/3+10/3=12','nonzero_divisible_by3':'9*(10/3)=30'},
    'center_density_formula':'3^17/sqrt(10)','center_density_squared':str(s.Rational(3**34,10)),
    'center_density_approx':float(s.Integer(3)**17/s.sqrt(10)),
    'world_record_status':'priority not established',
    'independent_CVP_certificate':certificate_info,
}
(P/'ne65_verified.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
