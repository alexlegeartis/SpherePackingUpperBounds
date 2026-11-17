from fractions import Fraction as F
from math import factorial
# Positive Taylor partial sums are rigorous lower bounds on exp(x).
def exp_lower(x,n):return sum((F(x)**k/F(factorial(k))for k in range(n+1)),F(0))
R_upper=F(1,2*10**18)
assert exp_lower(120,400)>(F(3,2)**192)/R_upper
q=4*10**8
compact=81*q*R_upper*(121+F(22,56)+F(2,56**2))
assert compact<F(2,10**6)
ani_poly=F(38,3)**2+F(38,36)+F(1,288)
assert exp_lower(64,250)>81*q*ani_poly*10**15
assert exp_lower(240,700)>q*3**96*10**45
print('PASS exact exponential/rational constants: compact tail <2e-6; anisotropic positive tail <1e-15; q*h(3)<1e-45.')
print('Requires geometric counting lemmas and completed compact main gate; NOT a standalone density certificate.')
