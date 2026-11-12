"""Exact shared count-to-radial gate; safe to import without running proof."""
from fractions import Fraction as F
def count_linkage(objective,radial_count):
 objective=F(objective);floor=objective.numerator//objective.denominator
 assert floor==7724
 assert F(radial_count)==floor
 return floor
