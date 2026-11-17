from pathlib import Path
import json,sympy as s,hashlib
from fractions import Fraction as F
p=Path(__file__).resolve().parent; c=json.load(open(p/'certificate.json'));assert c==json.load(open(p.parent/'combined48_trial/certificate.json'))
x=s.symbols('x');M=s.Matrix(c['fourier_sos']);N=s.Matrix(c['fourier_halfline_sos']);b=s.Matrix([(x-1)**i for i in range(5)]);d=s.Matrix([(x-1)**i for i in range(4)]);m=s.Poly((b.T*M*b)[0]+x*(d.T*N*d)[0],x)
assert all(m.nth(i)==s.Rational(v)for i,v in enumerate(c['fourier_additive_coefficients']))
assert all(M[:k,:k].det()>0 for k in range(1,6));assert all(N[:k,:k].det()>0 for k in range(1,5));assert F(c['density'])<=F('3.527782891e-7')
print('Direct symbolic SOS reconstruction, positive definiteness, full certificate equality and safe rounding PASS')
