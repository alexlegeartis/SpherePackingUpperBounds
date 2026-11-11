"""Independent rational Schur factorization and outward density rounding.

Run with an optional repository root argument. Without an argument, expects
installation at verifications/cyclic_origin_slack/verify_rational.py.
Source continuum checks remain dependencies of the package's main verifier.
"""
import json,math,sys
from fractions import Fraction as F
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions must be enabled')
root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[2]

def atan_bounds(x):
    assert 0<x<1
    lo=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(32))
    return lo,lo+x**65/65

# Machin identity pi=16 atan(1/5)-4 atan(1/239), with alternating remainders.
l5,u5=atan_bounds(F(1,5));l239,u239=atan_bounds(F(1,239))
pi_upper=16*u5-4*l239
for n,threshold in [(18,'0.013224473'),(19,'0.009521000'),(20,'0.006936042')]:
    d=json.loads((root/f'verifications/dim{n}/refinement_cyclic/certificate.json').read_text())
    assert d['dimension']==n
    A=[[F(x) for x in row] for row in d['matrices'][0]]
    size=len(A);assert all(len(row)==size for row in A)
    assert all(A[i][j]==A[j][i] for i in range(size) for j in range(size))
    m=size-1;B=[row[1:] for row in A[1:]];b=[row[0] for row in A[1:]]
    aug=[B[i][:]+[b[i]] for i in range(m)]
    # Positive elimination pivots certify B positive definite exactly.
    for k in range(m):
        pivot=aug[k][k];assert pivot>0
        for i in range(k+1,m):
            factor=aug[i][k]/pivot
            for j in range(k,m+1):aug[i][j]-=factor*aug[k][j]
    x=[F(0)]*m
    for i in reversed(range(m)):
        x[i]=(aug[i][m]-sum(aug[i][j]*x[j] for j in range(i+1,m)))/aug[i][i]
    assert all(sum(B[i][j]*x[j] for j in range(m))==b[i] for i in range(m))
    alpha=A[0][0]-sum(b[i]*x[i] for i in range(m));assert alpha>0
    W=[x]+[[F(i==j) for j in range(m)] for i in range(m)]
    assert all(A[i][j]-(alpha if i==j==0 else 0)==sum(W[i][k]*B[k][ell]*W[j][ell] for k in range(m) for ell in range(m)) for i in range(size) for j in range(size))
    norm=F(d['normalization']);assert norm>0
    base=(F(d['physical_polynomial'][0])+A[0][0])/norm
    assert base==F(d['exact_objective'])
    objective=base-alpha/norm;s=F(d['exclusion_squared']);assert objective>0 and s>0
    # Ball radius sqrt(s)/2. Square the density to eliminate odd-dimensional sqrt(s).
    if n%2==0:
        upper=(objective*pi_upper**(n//2)*(s/4)**(n//2)/math.factorial(n//2))**2
    else:
        upper=objective**2*pi_upper**(n-1)*s**n*4*math.factorial((n+1)//2)**2/math.factorial(n+1)**2
    assert upper<F(threshold)**2
    print(f'PASS {n}D: exact Schur PSD factorization and rational Machin density <{threshold}',flush=True)
print('PASS all three independent rational saturation checks')
