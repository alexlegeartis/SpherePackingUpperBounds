if not __debug__:
 raise RuntimeError('Do not run with python -O')
from fractions import Fraction as F
S=F(4,15)
assert S<F(8,29)
assert F(41)*(1+S/32)**8<44
assert (1+F(5,4)*S)*F(3,4)==1
assert (1+F(5,4)*S)*F(8,11)<1
allowed=[i for i in range(33,44) if all(i%p for p in [5,7,11,13,17,19,23,29,31,37,41,43]) and i%9]
assert allowed==[]
print('Strong near-frame rational constants PASS')
