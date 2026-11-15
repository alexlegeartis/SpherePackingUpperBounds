"""Finite rational-coordinate raw V7 energy, with outward rational intervals.

An evaluator, not a search or a universal floor. Run from any folder; inputs come from paths.py.
The original broad capacity is 1358. Stored count repairs enter exactly once.
The complete unlifted pair is reconstructed from singleton energy minus A.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
import argparse
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import paths
CHECKS = HERE
SOURCE = paths.V7
POINTED = HERE
sys.set_int_max_str_digits(0)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


p = load("d411_pointed_helpers", CHECKS / "pointed_graph_energy.py")
g = load("d411_gram_helpers", CHECKS / "energy_d4_gram.py")
allocation = load("d411_allocation_helpers", POINTED / "two_capacity_antipodal_allocation_exact.py")
pt, add, neg, mul, div, scale, isum = p.pt, p.add, p.neg, p.mul, p.div, p.scale, p.isum
kernel = json.loads((SOURCE / "kernel.json").read_text())
origin = json.loads((SOURCE / "origin.json").read_text())
scalar = json.loads((SOURCE / "scalar.json").read_text())
graphs = json.loads((SOURCE / "graph_energies.json").read_text())
acceptance = json.loads(paths.ACCEPTANCE.read_text())
R2 = Q(49, 16)
DEN = int(kernel["common_denominator"])
RADII5 = tuple(map(Q, origin["conditional_count_polynomials"][0]["radii"]))
RADII3 = tuple(map(Q, origin["nested_count_quadratic"]["radii"]))
U_RAW = Q(acceptance["exact_source_objective"])
A = 8 * U_RAW - Q(scalar["F"][0])
assert kernel["dimension"] == 4 and kernel["angular_degree"] == 18
assert kernel["local_degree"] == 6 and len(origin["u"]) == 49
assert list(map(Q, graphs["weights"][:2])) == [0, 0]
assert all(Q(row["rounding_origin_repair"]) == 0 for row in origin["conditional_count_polynomials"])


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


@lru_cache(None)
def features(m):
    """49 Bernstein features; closed-inner convention at every break."""
    assert 1 <= m <= R2
    return tuple(p.features(m, kernel))


@lru_cache(None)
def amplitude(m, degree):
    return tuple(isum(scale(z, Q(v)) for z, v in zip(features(m), ray["vector"][int(degree == 0):]))
                 for ray in kernel["blocks"][degree]["rays"])


@lru_cache(None)
def zonal(m, n, h, degree):
    """U_degree(cos)/(degree+1), h=2 dot; Gegenbauer parameter 1."""
    assert h * h <= 4 * m * n
    t2 = h * h / (4 * m * n)
    if degree % 2 == 0:
        return pt(g.even_zonal(t2, degree))
    coefficient = sum((-1) ** j * comb(degree-j, j) * 2 ** (degree-2*j)
                      * t2 ** ((degree-2*j-1)//2) for j in range((degree+1)//2)) / Q(degree+1)
    return scale(div(pt(h), scale(p.root(m*n), 2)), coefficient)


def count_energy(counts3, counts5):
    """The stored constants ALREADY include their rounding repairs."""
    value = g.polynomial(origin["nested_count_quadratic"]["coefficients"], counts3)
    return value + sum(Q(row["weight"]) * g.polynomial(row["coefficients"], counts5)
                       for row in origin["conditional_count_polynomials"])


def radial_weight(m, row, cutoff=True):
    if cutoff and "pivot_lower_cutoff" in row and m <= Q(row["pivot_lower_cutoff"]) ** 2:
        return pt(0)
    return isum(scale(z, Q(c)) for z, c in zip(features(m), row["w"]))


def graph_energy(hist):
    specs = [
        ("broad_antipodal", Q(graphs["weights"][2]), 1358, Q(33,28), Q(7,4), Q(-1), Q(-7,8)),
        ("near_antipodal", Q(graphs["outer_graph_weight"]), 120, Q(1), Q(103,100), Q(-1), Q(-3,4)),
        ("narrow_antipodal", Q(graphs["narrow_graph_weight"]), 604, Q(33,28), Q(7,4), Q(-1), Q(-183,200)),
        ("moderate", Q(graphs["moderate_graph_weight"]), 1880, Q(7,5), Q(3,2), Q(-4,5), Q(-11,20)),
    ]
    rows = []
    for name, weight, capacity, low_r, high_r, low_t, high_t in specs:
        edges = sum(count for (m,n,h),count in hist.items()
                    if low_r**2 <= m <= high_r**2 and low_r**2 <= n <= high_r**2
                    and h < 0 and 4*m*n*low_t**2 >= h*h >= 4*m*n*high_t**2)
        assert 0 <= edges <= capacity
        rows.append(dict(name=name, edges=edges, capacity=capacity, energy=weight*(capacity-edges)))
    return sum(row["energy"] for row in rows), rows


def unlifted_components(shells, hist):
    """Energy components before the full-pair count gauge.

    shells maps squared radius to multiplicity; hist maps (m,n,2 dot) to
    ordered multiplicity, including diagonal pairs. Coordinates are not
    needed after exact separation and incidence have been established.
    """
    assert all(1 <= m <= R2 and isinstance(count, int) and count > 0 for m,count in shells.items())
    assert sum(hist.values()) == sum(shells.values()) ** 2
    z = [isum(scale(features(m)[i], count) for m,count in shells.items()) for i in range(49)]
    e0 = pt(0)
    for ray in kernel["blocks"][0]["rays"]:
        assert ray["vector"][0] == 0 and ray["weight"] > 0
        value = isum(scale(add(zz, pt(Q(u))), Q(v)) for zz,u,v in zip(z,origin["u"],ray["vector"][1:]))
        e0 = add(e0, scale(g.square(value), Q(ray["weight"], DEN)))
    angular = {}
    for degree,block in enumerate(kernel["blocks"][1:],1):
        harmonic = {(m,n): isum(scale(zonal(m,n,h,degree),count)
                    for (mm,nn,h),count in hist.items() if (mm,nn)==(m,n)) for m in shells for n in shells}
        value = pt(0)
        for i,ray in enumerate(block["rays"]):
            term = isum(mul(mul(amplitude(m,degree)[i],amplitude(n,degree)[i]),harmonic[m,n])
                        for m in shells for n in shells)
            value = add(value,scale(term,Q(ray["weight"],DEN)))
        # The inherited PSD Gram-ray theorem supplies exact nonnegativity.
        assert value[1] >= 0
        angular[degree] = (max(Q(0),value[0]),value[1])
    counts5 = [sum(count for m,count in shells.items() if m <= r*r) for r in RADII5]
    counts3 = [sum(count for m,count in shells.items() if m <= r*r) for r in RADII3]
    spherical = pt(0)
    for row in origin["appended_spherical_local_energy_rays"]:
        limit = Q(row["R"])**2
        count = sum(count for m,count in shells.items() if m <= limit)
        formed = isum(scale(isum(scale(zonal(m,n,h,degree),Q(c))
                      for degree,c in enumerate(row["scalar_coefficients"])),mult)
                      for (m,n,h),mult in hist.items() if m <= limit and n <= limit)
        spherical = add(spherical,scale(add(pt(count*Q(row["B"])),neg(formed)),Q(row["weight"])))
    caps = pt(0)
    for row in origin["pointed_caps"]:
        for n,count in shells.items():
            occupied = isum(scale(p.probability(m,n,h,row),mult)
                            for (m,nn,h),mult in hist.items() if nn==n)
            gap = add(pt(Q(row["M"])*count),neg(occupied))
            caps = add(caps,scale(mul(radial_weight(n,row),gap),Q(row["weight"])))
    tents = pt(0)
    for row in origin["pointed_tents"]:
        d = Q(row["d"])
        assert Q(row["sigma"]) == -1
        for n,count in shells.items():
            masses = []
            for (m,nn,h),mult in hist.items():
                if nn != n:
                    continue
                raw = add(pt(1-2*(m+d*d)),neg(div(pt(2*d*h),p.root(n))))
                masses.append(scale((max(Q(0),raw[0]),max(Q(0),raw[1])),mult))
            gap = add(pt(Q(row["capacity"])*count),neg(isum(masses)))
            tents = add(tents,scale(mul(radial_weight(n,row),gap),Q(row["weight"])))
    allocations = pt(0)
    for row in origin["independent_antipodal_allocations"]:
        for n,count in shells.items():
            occupied = isum(scale(allocation.capacity_point_bounds(m+n+h,row["capacity_profile"],bits=160),mult)
                            for (m,nn,h),mult in hist.items() if nn==n)
            gap = add(pt(count),neg(occupied))
            allocations = add(allocations,scale(mul(radial_weight(n,row),gap),Q(row["weight"])))
    conditional = Q(0)
    for row in origin["conditional_radial_rays"]:
        ni = sum(count for m,count in shells.items() if m <= Q(row["inner_radius"])**2)
        no = sum(count for m,count in shells.items() if m <= Q(row["outer_radius"])**2)
        conditional += Q(row["weight"])*(23-ni)*(no-ni)
    graph, graph_rows = graph_energy(hist)
    components = dict(E0=e0,angular=isum(angular.values()),count=pt(count_energy(counts3,counts5)),
                      spherical=spherical,caps=caps,tents=tents,allocations=allocations,
                      conditional=pt(conditional),graphs=pt(graph))
    return components,dict(counts3=counts3,counts5=counts5,angular_by_degree=angular,graphs=graph_rows)


@lru_cache(None)
def scalar_value(m):
    x = Q(scalar["alpha"])*m
    assert 0 <= x <= 7
    term = total = Q(1)
    for k in range(1,221):
        term *= x/k
        total += term
    tail = term*x/221/(1-x/222)
    polynomial = Q(0)
    for c in reversed(scalar["F"]):
        polynomial = polynomial*x+Q(c)
    return scale((1/(total+tail),1/total),polynomial)


@lru_cache(None)
def unlifted_pair(m):
    """q0=f+2h+d from the exact singleton identity E({y})-E(empty)."""
    components,_ = unlifted_components({m:1},Counter({(m,m,2*m):1}))
    result = add(scalar_value(m),add(isum(components.values()),pt(-A)))
    assert result[0] <= 0  # fail closed on a contradiction to accepted signs
    return result[0],min(Q(0),result[1])


def evaluate(points, anchor=None):
    """Return complete raw energy for one anchor in a finite rational set.

    All input centres are checked for unit separation. Points outside the
    closed support ball are retained for separation but have no energy term.
    Passing a truncated set evaluates that finite packing, never a completion.
    """
    points = [tuple(map(Q,row)) for row in points]
    assert all(len(row)==4 for row in points) and len(set(points))==len(points)
    anchor = tuple(map(Q,anchor)) if anchor is not None else (Q(0),)*4
    assert anchor in points
    assert all(dot(tuple(a-b for a,b in zip(x,y)),tuple(a-b for a,b in zip(x,y))) >= 1
               for x,y in combinations(points,2))
    vectors = [tuple(x-a for x,a in zip(row,anchor)) for row in points if row!=anchor]
    active = [v for v in vectors if dot(v,v) <= R2]
    shells = Counter(dot(v,v) for v in active)
    hist = Counter((dot(x,x),dot(y,y),2*dot(x,y)) for x in active for y in active)
    components,details = unlifted_components(shells,hist)
    gauge_count = sum(count for m,count in shells.items() if m <= Q(origin["full_pair_count_gauge"]["radius"])**2)
    assert gauge_count <= 24
    gauge = pt(0)
    if gauge_count < 24:
        weights = []
        for m,count in shells.items():
            psi = Q(int(m<=Q(origin["full_pair_count_gauge"]["radius"])**2),24)
            weights.append(scale(neg(unlifted_pair(m)),Q(count)/(1-psi)))
        gauge = scale(isum(weights),1-Q(gauge_count,24))
    components["gauge"] = gauge
    total = isum(components.values())
    assert total[0] <= total[1]
    return dict(status="EXACT_FINITE_RAW_V7_ENCLOSURE",dimension=4,gegenbauer_parameter="1",
                raw_convention="original broad capacity 1358; stored count repair once",
                anchor=anchor,point_count=len(points),support_count=len(active),shell_counts=dict(shells),
                gauge_count=gauge_count,origin_constant=A,components=components,total=total,**details)


def encode(value):
    if isinstance(value,Q):
        return str(value)
    if isinstance(value,dict):
        return {str(key):encode(item) for key,item in value.items()}
    if isinstance(value,tuple) and len(value)==2 and all(isinstance(item,Q) for item in value):
        # Compact outward serialization; calculations above keep full rationals.
        denominator = 1 << 128
        low,high = value
        lower = Q(low.numerator*denominator//low.denominator,denominator)
        upper = Q(-((-high.numerator*denominator)//high.denominator),denominator)
        assert lower <= low <= high <= upper
        return [str(lower),str(upper)]
    if isinstance(value,(tuple,list)):
        return [encode(item) for item in value]
    return value


def identity_check():
    components,_ = unlifted_components({},Counter())
    assert components["E0"] == pt(Q(origin["kappa"]))
    assert isum(components.values()) == pt(A)
    assert A == Q(origin["kappa"])+count_energy([0]*3,[0]*5)+Q(4769,25000000)
    return dict(status="EXACT_RAW_V7_ORIGIN_IDENTITY_PASS",origin_constant=A,
                components=components,count_repair_already_serialized=Q(origin["nested_count_quadratic"]["constant_rounding_repair"]),
                raw_objective=U_RAW,target_credit=8*(U_RAW-Q(629,1000)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,help="JSON object: points, optional anchor; rational strings or integers")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--identity-only",action="store_true")
    args = parser.parse_args()
    assert bool(args.input) != args.identity_only
    if args.identity_only:
        result = identity_check()
    else:
        data = json.loads(args.input.read_text())
        result = evaluate(data["points"],data.get("anchor"))
    args.output.write_text(json.dumps(encode(result),indent=2)+"\n")
    print(json.dumps({"status":result["status"],"output":str(args.output)}))
