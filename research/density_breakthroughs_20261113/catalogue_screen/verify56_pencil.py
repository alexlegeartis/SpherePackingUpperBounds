#!/usr/bin/env python3
"""Exact certificate excluding a new 56D record in the entire discovered pencil.

The elementary one-variable argument is given in README.md. This verifier
checks the matrix identities, integral equivalences, two short witnesses, and
the exact comparison with the current 56D benchmark. It does not use SVP.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
import sympy as s

P=Path(__file__).resolve().parent
cert=P/'L56a_pencil_certificate.npz'
if cert.exists():
    a=np.load(cert)
    G=s.Matrix(a['G']);D=s.Matrix(a['D']);U=s.Matrix(a['U']);R=s.Matrix(a['R'])
    v=s.Matrix(a['v']);w=s.Matrix(a['w'])
else:
    a=np.load(P/'L56a_flex.npz');o=np.load(P/'L56a_orbit.npz')
    G=s.Matrix(a['G']);D=s.Matrix(a['D']);B=29*G.inv()*D;U=(5*s.eye(56)+B)/2
    R=s.Matrix(o['generators'][2]).T
    v=s.eye(56)[:,0]
    witness=json.loads((P/'test56_3/report.json').read_text())['witnesses'][0]['vector']
    w=s.Matrix(witness)
    assert all(x.q==1 for x in U)
    np.savez_compressed(cert,**{k:np.array(x,dtype=np.int64)for k,x in {'G':G,'D':D,'U':U,'R':R,'v':v,'w':w}.items()})
n=56;I=s.eye(n)
assert G==G.T and D==D.T
assert G.det()==29**28
L,diag=G.LDLdecomposition(hermitian=False)
assert all(diag[i,i]>0 for i in range(n))
assert L*diag*L.T==G
H=G.inv()*D
assert H*H==I/29 and s.trace(H)==0
assert U.det()==1 and R.det() in (-1,1)
assert U.T*G*U==(27*G+145*D)/2
assert U.T*D*U==(5*G+27*D)/2
assert R.T*G*R==G and R.T*D*R==-D
assert (v.T*G*v)[0]==28 and (v.T*D*v)[0]==0
assert (w.T*G*w)[0]==32 and (w.T*D*w)[0]==-4
# On 0<=t<=1 the bound is [49/(29-t^2)]^14, maximized at t=1.
# On 1<=t<=5 it is [(8-t)^2/(29-t^2)]^14; derivative changes
# from negative to positive at29/8, so its maximum is at an endpoint.
assert F(1)<F(29,8)<F(5)
bound=max(F(7,4),F(9,4))**14
record=F(5,2)**28/F(29)**4
assert bound<record
report={'dimension':56,'determinant':'(29-t^2)^28','positive_definite_domain':'|t|<sqrt(29)','fundamental_parameter_interval':'0<=t<=5, by an integral Pell unit and reflection','minimum_upper_witnesses':['28','32-4t'],'entire_pencil_center_density_upper':str(bound),'entire_pencil_center_density_upper_decimal':float(bound),'record_center_density':str(record),'record_center_density_decimal':float(record),'ratio_to_record':str(bound/record),'ratio_to_record_decimal':float(bound/record),'new_record_in_pencil':False,'verification':'exact rational and integer matrix arithmetic; global in this one-parameter family'}
(P/'L56a_pencil_verified.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
