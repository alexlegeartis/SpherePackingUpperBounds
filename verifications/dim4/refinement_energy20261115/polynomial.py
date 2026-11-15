"""Exact arithmetic for normalized S^3 Gegenbauer energy polynomials."""
from fractions import Fraction as Q
from math import comb


def add(a, b):
    out = [Q(0)] * max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    return out


def power_coefficients(zonal):
    """Convert sum c_d U_d(t)/(d+1) to exact monomial coefficients."""
    u = [[Q(1)], [Q(0), Q(2)]]
    for d in range(2, len(zonal)):
        u.append(add([Q(0)] + [2*x for x in u[-1]], [-x for x in u[-2]]))
    out = [Q(0)] * len(zonal)
    for d, c in enumerate(zonal):
        for k, v in enumerate(u[d]): out[k] += c*v/(d+1)
    return out


def bernstein_on(power, lo, hi):
    """Exact coefficients on t=lo+(hi-lo)x, 0<=x<=1."""
    n = len(power)-1
    shifted = [sum(power[k]*comb(k,j)*lo**(k-j)*(hi-lo)**j
                   for k in range(j,n+1)) for j in range(n+1)]
    return [sum(shifted[j]*Q(comb(i,j),comb(n,j)) for j in range(i+1))
            for i in range(n+1)]


def check_nonpositive(zonal, a, cap):
    """Prove p<=0 throughout [-a,cap] by exact Bernstein subdivision."""
    neg = [-c for c in power_coefficients(zonal)]
    pending = [(-a,Q(0),0),(Q(0),cap,0)]
    leaves = 0
    minimum = None
    while pending:
        lo,hi,depth = pending.pop()
        bc = bernstein_on(neg,lo,hi)
        if min(bc)>=0:
            leaves += 1
            minimum=min(bc) if minimum is None else min(minimum,min(bc))
        else:
            if depth>=60:
                raise ValueError('central sign unresolved: '+str((lo,hi)))
            mid=(lo+hi)/2
            pending.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
    return dict(leaves=leaves,minimum_bernstein_coefficient=str(minimum),
                interval=[str(-a),str(cap)])


def energy(zonal, n):
    return n*n*zonal[0]-n*sum(zonal)
