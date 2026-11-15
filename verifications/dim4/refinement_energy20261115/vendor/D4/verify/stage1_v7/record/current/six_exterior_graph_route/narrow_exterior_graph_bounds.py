"""Exact conservative subtraction for the separately proved narrow degree4 antipodal graph."""
from fractions import Fraction as Q
from graph_repair_bounds import contained,cosine_exterior

CAPACITY=604
DEFAULT_WEIGHT=Q(3,12500000)

def cyclic_decrement_lower(kind,bounds,weight=DEFAULT_WEIGHT):
    assert isinstance(weight,Q) and weight>=0
    if kind=='full':return Q(0)
    assert kind=='outside' and len(bounds)==3
    a,b,z=[tuple(map(Q,box)) for box in bounds]
    outer=(Q(33,28),Q(7,4));angular=(Q(-1),Q(-183,200))
    return weight*int(contained(a,outer) and contained(b,outer) and contained(cosine_exterior(a,b,z),angular))

def density_cost(weight=DEFAULT_WEIGHT):
    assert isinstance(weight,Q) and weight>=0
    return CAPACITY*weight/8
