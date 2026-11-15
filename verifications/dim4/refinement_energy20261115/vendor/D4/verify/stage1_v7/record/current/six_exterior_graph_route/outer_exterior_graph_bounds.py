"""Exact closed-box lower/upper corrections for the near antipodal graph.

Local graph positivity and its origin charge require a separate fresh
verify_near_antipodal_exterior_graph.py proof.  This module only evaluates
the indicator kernel and its conservative box enclosure.
"""
from fractions import Fraction as Q

CAPACITY = 120
RADIAL = (Q(1),Q(103,100))
ANGULAR = (Q(-1),Q(-3,4))
KERNEL_SUPPORT_RADIUS = Q(7,4)


def _cosine_bounds(bounds):
    assert len(bounds)==3
    (rlo,rhi),(slo,shi),(zlo,zhi)=bounds
    assert 1<=rlo<=rhi<=KERNEL_SUPPORT_RADIUS
    assert 1<=slo<=shi<=KERNEL_SUPPORT_RADIUS and 0<=zlo<=zhi<=1
    def cosine(r,s,z):
        return -1+z*((r+s)**2-KERNEL_SUPPORT_RADIUS**2)/(2*r*s)
    # The expression increases in z. Its r derivative has the sign of
    # z*(r^2-s^2+R^2), which is nonnegative on 1<=r,s<=R; likewise in s.
    answer=cosine(rlo,slo,zlo),cosine(rhi,shi,zhi)
    assert -1<=answer[0]<=answer[1]<=1
    return answer


def cyclic_decrement_lower(kind,bounds,weight):
    assert kind in ('full','outside') and type(weight) is Q and weight>=0
    if kind=='full' or weight==0:
        return Q(0)
    assert len(bounds)==3
    cosine=_cosine_bounds(bounds)
    if all(RADIAL[0]<=lo<=hi<=RADIAL[1] for lo,hi in bounds[:2]) and cosine[1]<=ANGULAR[1]:
        return weight
    return Q(0)


def cyclic_decrement_upper(kind,bounds,weight):
    assert kind in ('full','outside') and type(weight) is Q and weight>=0
    if kind=='full' or weight==0:
        return Q(0)
    assert len(bounds)==3
    cosine=_cosine_bounds(bounds)
    if any(hi<RADIAL[0] or lo>RADIAL[1] for lo,hi in bounds[:2]) or cosine[0]>ANGULAR[1]:
        return Q(0)
    return weight


def density_cost(weight):
    assert type(weight) is Q and weight>=0
    return CAPACITY*weight/8
