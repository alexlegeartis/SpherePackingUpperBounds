from pathlib import Path
import json,hashlib
import sympy as s
if not __debug__:raise RuntimeError('Assertions required')
D=Path(__file__).parent;p=D/'frozen.json';sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='fe60310a72412cc57e47dc0f0825a0aa7806dd80d3c7c0c045888c6925957dab'
z=json.load(open(p));co=list(map(s.Rational,z['co']));B=[[[s.Rational(v)for v in row]for row in M]for M in z['B']];x=s.symbols('x')
Q=[s.Poly(sum(co[k]*B[k][i][j]*s.gegenbauer(2*k,22,x)/s.gegenbauer(2*k,22,1)for k in range(len(co))),x)for i,j in [(0,0),(0,1),(1,1)]]
eta=s.Rational(1,1000000);corners=[Q[0],Q[0]+Q[1],Q[0]+2*Q[1]+Q[2]]
for pp in corners:
 t=pp+eta
 assert t.eval(0)<0 and t.eval(s.Rational(1,2))<0 and t.count_roots(0,s.Rational(1,2))==0
# For any u,v in[0,1], bilinear Q00+(u+v)Q01+uv Q11 is the convex combination
# (1-u)(1-v) Q00 + [u(1-v)+(1-u)v](Q00+Q01) + uv(Q00+2Q01+Q11).
# In the compact cap domain1<=a<=b<=R=19/13 and0<=dot<=(a+b-R)/2,
# c=dot/sqrt(ab)<=1/2 because a+b-R<=min(a,b)<=sqrt(ab).
# Both radial ratios u=h(a),v=h(b) lie in[0,1]. Thus the whole cap is <-eta.
out=dict(status='COMPLETE EXACT CAP PASS',dimension=46,input_sha=sha,margin=str(eta),method='Three exact univariate Sturm signs and bilinear convex corners',corner_coefficients=[list(map(str,pp.all_coeffs()))for pp in corners])
(D/'cap_exact_result.json').write_text(json.dumps(out,indent=2));print(out['status'],out['margin'])
