import sys,json
sys.path.insert(0,'solver_packages')
from flint import fmpq,fmpq_poly
from fractions import Fraction as F
j=json.load(open('continuum_three_point/cyclic_angular40_discovery_snapshot.json'));cs=[fmpq(F(str(v)).numerator,F(str(v)).denominator)for v in j['coefficients']];eps=fmpq(1,10000);cs[0]=1-eps;cs=[v/(1-eps)for v in cs];assert min(cs)>=0
X=fmpq_poly([0,1]);G=[fmpq_poly([1]),X]
for k in range(1,len(cs)-1):G.append(((2*k+46)*X*G[-1]-k*G[-2])/(k+46))
Q=sum((a*g for a,g in zip(cs,G)),fmpq_poly([]))
def roots(P,a,b):
 s=[P,P.derivative()]
 while s[-1].degree()>0:
  z=(-(s[-2]%s[-1])).numer();s.append(fmpq_poly(z/abs(z.content())))
 def V(t):
  sg=[1 if x>0 else -1 for x in [p(t)for p in s]if x];return sum(a!=b for a,b in zip(sg,sg[1:]))
 return V(a),V(b)
c=fmpq(5001,10000);cap=Q+fmpq(1,100000);v=roots(cap,-1,fmpq(1,2));assert v[0]==v[1]and cap(-1)<0 and cap(fmpq(1,2))<0
v0=roots(Q,-1,c);v1=roots(Q,c,1);assert v0[0]==v0[1]and v1[0]-v1[1]==1 and Q(c)<0 and Q(1)>0
out={'dimension':48,'degree':Q.degree(),'alpha':40,'cutoff':str(c),'gegenbauer_coefficients':list(map(str,cs)),'polynomial_coefficients':list(map(str,Q.coeffs())),'B_exact':str(sum(cs)),'B_float':float(sum(cs)),'subtracted_epsilon':str(eps),'negative_cap_delta':'1/100000','cap_variations':v,'single_positive_root_variations':v1}
json.dump(out,open('continuum_three_point/cyclic_angular40_exact.json','w'),indent=2);print('PASS angular',out['B_float'])
