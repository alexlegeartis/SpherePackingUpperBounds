"""Exact lower and conservative upper for the moderate single-color graph."""
from fractions import Fraction as Q
from graph_repair_bounds import contained,cosine_exterior

CAPACITY=1880
RADIAL=(Q(7,5),Q(3,2))
ANGULAR=(Q(-4,5),Q(-11,20))
KERNEL_SUPPORT_RADIUS=Q(7,4)


def cyclic_decrement_lower(kind,bounds,weight=Q(0)):
    assert type(weight) is Q and weight>=0
    assert kind in ('full','outside') and len(bounds)==3
    if kind=='full':return Q(0)
    a,b,z=[tuple(map(Q,interval)) for interval in bounds]
    return weight*int(contained(a,RADIAL) and contained(b,RADIAL)
        and contained(cosine_exterior(a,b,z),ANGULAR))


def cyclic_decrement_upper(kind,bounds,weight=Q(0)):
    assert type(weight) is Q and weight>=0
    assert kind in ('full','outside') and len(bounds)==3
    if kind=='full':return Q(0)
    if all(Q(left)==Q(right) for left,right in bounds):
        return cyclic_decrement_lower(kind,bounds,weight)
    return weight


def density_cost(weight=Q(0)):
    assert type(weight) is Q and weight>=0
    return CAPACITY*weight/8
