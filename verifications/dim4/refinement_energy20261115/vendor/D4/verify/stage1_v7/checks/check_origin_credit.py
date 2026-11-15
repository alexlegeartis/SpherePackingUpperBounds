#!/usr/bin/env python3
"""Independent exact audit of the retained D4 V7 radial origin credit.

Uses only the Python standard library. Reads the original files without writing
to them. It reconstructs the matrix from integer Gram rays, verifies all vectors,
proves all seven radial signs by rational de Casteljau subdivision, and compares
the fully charged rational objective to the requested strict threshold.

This checks the origin-credit and objective subcertificate. It does not prove
the scalar, graph capacities, pair, or cyclic global prerequisites.
"""
from fractions import Fraction as Q
from pathlib import Path
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math


def require(test, reason):
    if not test:
        raise ValueError(reason)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact(value):
    require(type(value) in (str, int), "An exact scalar must be a string or integer")
    return Q(value)


def dot(a, b):
    require(len(a) == len(b), "Vector dimension mismatch")
    return sum((x*y for x, y in zip(a, b)), Q(0))


def split(coefficients, t):
    """Exact de Casteljau subdivision at t in [0,1]."""
    require(Q(0) <= t <= Q(1), "Invalid subdivision parameter")
    rows = [list(coefficients)]
    while len(rows[-1]) > 1:
        old = rows[-1]
        rows.append([(1-t)*x+t*y for x, y in zip(old, old[1:])])
    return [row[0] for row in rows], [row[-1] for row in reversed(rows)]


def restrict(coefficients, left, right):
    require(Q(0) <= left < right <= Q(1), "Invalid closed interval")
    suffix = split(coefficients, left)[1]
    return split(suffix, (right-left)/(1-left))[0]


def power_from_bernstein(coefficients):
    degree = len(coefficients)-1
    return [Q(math.comb(degree, k))*sum(
        (Q((-1)**(k-j)*math.comb(k, j))*coefficients[j]
         for j in range(k+1)), Q(0)) for k in range(degree+1)]


def sign_certificate(coefficients, left, right):
    queue = [(left, right, restrict(coefficients, left, right), 0)]
    leaves = []
    nodes = 0
    while queue:
        a, b, values, depth = queue.pop()
        nodes += 1
        require(nodes <= 100000, "Radial sign subdivision budget exceeded")
        if min(values) >= 0:
            leaves.append({"left": str(a), "right": str(b),
                           "minimum_bernstein_coefficient": str(min(values)),
                           "bernstein_coefficients": list(map(str, values)),
                           "depth": depth})
            continue
        require(max(values) >= 0, "Radial weight is negative on a subinterval")
        require(depth < 128, "Radial sign depth budget exceeded")
        lc, rc = split(values, Q(1, 2))
        mid = (a+b)/2
        queue.append((mid, b, rc, depth+1))
        queue.append((a, mid, lc, depth+1))
    leaves.sort(key=lambda leaf: Q(leaf["left"]))
    require(Q(leaves[0]["left"]) == left and Q(leaves[-1]["right"]) == right,
            "Radial proof omits an endpoint")
    require(all(x["right"] == y["left"] for x, y in zip(leaves, leaves[1:])),
            "Radial proof leaves a coverage gap")
    return {"method": "EXACT_RATIONAL_BERNSTEIN_SUBDIVISION",
            "closed_interval_covered": True, "nodes": nodes,
            "terminal_boxes": len(leaves), "leaves": leaves}


def dec(value):
    with localcontext() as ctx:
        ctx.prec = 70
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def audit(archive, output):
    archive = archive.resolve()
    output = output.resolve()
    require(archive.exists(), "Archive directory is missing")
    require(output != archive and not output.is_relative_to(archive),
            "Audit output must be outside the immutable original archive")
    candidate = "current/candidates/round9_six_graphs_degree18"
    paths = {
        "kernel": candidate+"/kernel.json",
        "origin": candidate+"/origin.json",
        "scalar": candidate+"/scalar.json",
        "graphs": candidate+"/graph_energies.json",
        "credit": "audit/round9_origin_credit/ORIGIN_CREDIT_VERIFIED.json",
        "composition": "audit/round9_six_graphs_monotone_origin_credit_global_verified.json",
    }
    hashes = {p: sha(archive/p) for p in paths.values()}
    data = {key: json.loads((archive/path).read_text()) for key, path in paths.items()}
    kernel, origin, scalar = data["kernel"], data["origin"], data["scalar"]
    credit, graphs, composition = data["credit"], data["graphs"], data["composition"]
    require(set(credit["source_sha256"]) == {"kernel.json", "origin.json", "scalar.json"},
            "Unexpected credit source binding")
    for name, digest in credit["source_sha256"].items():
        require(sha(archive/candidate/name) == digest,
                "Retained origin credit is not bound to the V7 "+name)
    require(graphs["base_input_sha256"] == credit["source_sha256"],
            "Graph model and origin credit use different base files")
    require(exact(scalar["alpha"]) == 2 and exact(scalar["G"][0]) == 1,
            "Unexpected density normalization")
    degree = kernel["local_degree"]
    breaks = list(map(exact, kernel["breaks"]))
    intervals = [list(map(exact, pair)) for pair in kernel["radial_basis_intervals"]]
    nc = degree+1
    n = nc*(len(breaks)-1)
    require(degree == 6 and n == 49 and len(intervals) == 7,
            "Unexpected radial feature model")
    require(all(x < y for x, y in zip(breaks, breaks[1:])), "Unordered radial breaks")
    denominator = int(kernel["common_denominator"])
    require(denominator > 0, "Nonpositive Gram denominator")
    block = kernel["blocks"][0]
    require(block["dimension"] == n+1, "Unexpected degree-zero Gram dimension")
    mint = [[0]*n for _ in range(n)]
    for ray in block["rays"]:
        w, vector = ray["weight"], ray["vector"]
        require(type(w) is int and w > 0 and len(vector) == n+1,
                "Invalid positive integer Gram ray")
        require(all(type(v) is int for v in vector), "Noninteger Gram ray")
        for i, x in enumerate(vector[1:]):
            for j, y in enumerate(vector[1:]):
                mint[i][j] += w*x*y
    require(all(mint[i][j] == mint[j][i] for i in range(n) for j in range(n)),
            "Nonsymmetric reconstructed Gram matrix")
    def mat(vector):
        return [dot(row, vector)/denominator for row in mint]
    u, v = (list(map(exact, origin[key])) for key in ("u", "v"))
    a = list(map(exact, credit["weight_Bernstein_coefficients"]))
    b = list(map(exact, credit["range_solution"]))
    require(len(u) == len(v) == len(a) == len(b) == n, "Invalid origin vector dimension")
    require(mat(u) == v, "M u differs from v")
    require(mat(b) == a, "M b differs from a")
    kappa = exact(origin["kappa"])
    require(dot(u, v) == kappa, "Original completion constant is inconsistent")
    gain = 2*dot(u, a)-dot(a, b)
    require(gain > 0, "Nonpositive origin credit")
    unew = [x-y for x, y in zip(u, b)]
    vnew = [x-y for x, y in zip(v, a)]
    require(mat(unew) == vnew, "Shifted range equations fail")
    require(dot(unew, vnew) == kappa-gain >= 0, "Shifted PSD origin constant fails")
    require(unew == list(map(exact, credit["new_psd_origin_u"])), "Retained shifted u mismatch")
    require(vnew == list(map(exact, credit["new_psd_origin_v"])), "Retained shifted v mismatch")
    require(exact(credit["new_psd_kappa"]) == kappa-gain, "Retained shifted kappa mismatch")
    require(exact(credit["origin_credit"]) == gain and exact(credit["density_credit"]) == gain/8,
            "Retained gain mismatch")
    signs = []
    for region, (lo, hi) in enumerate(intervals):
        left, right = breaks[region:region+2]
        require(lo <= left < right <= hi, "Physical interval lies outside mother interval")
        local = a[nc*region:nc*(region+1)]
        l, r = (left-lo)/(hi-lo), (right-lo)/(hi-lo)
        proof = sign_certificate(local, l, r)
        proof.update(region=region, physical_interval=list(map(str, (left, right))),
                     mother_interval=list(map(str, (lo, hi))),
                     normalized_interval=list(map(str, (l, r))),
                     original_bernstein_coefficients=list(map(str, local)),
                     original_power_coefficients=list(map(str, power_from_bernstein(local))))
        recorded = credit["radial_nonnegativity"][region]
        require(list(map(exact, recorded["power_coefficients"])) == power_from_bernstein(local),
                "Retained radial power polynomial does not match actual weight")
        signs.append(proof)
    capacities = [862, 1102, 1358, 120, 604, 1880]
    names = ["A", "B", "E9", "near_antipodal", "E4", "moderate"]
    require(graphs["schema"] == "INDEPENDENT_EXTERIOR_GRAPH_ENERGIES_V7", "Unexpected graph schema")
    weights = list(map(exact, graphs["weights"])) + [exact(graphs[key]) for key in
        ("outer_graph_weight", "narrow_graph_weight", "moderate_graph_weight")]
    require(len(weights) == len(capacities) and min(weights) >= 0, "Invalid graph weights")
    graph_charge = dot(capacities, weights)/8
    require(graph_charge == Q(4769, 200000000), "Unexpected exact total graph charge")
    base = exact(origin["objective"])
    raw = base+graph_charge
    objective = raw-gain/8
    threshold = Q(78749932299, 125000000000)
    require(exact(credit["source_objective"]) == base, "Credit base objective mismatch")
    require(exact(credit["new_formal_objective"]) == base-gain/8,
            "Credit base-only objective mismatch")
    require(exact(composition["exact_source_objective"]) == raw, "Composition raw objective mismatch")
    require(exact(composition["objective"]) == objective, "Composition final rational objective mismatch")
    require(exact(composition["density_credit"]) == gain/8, "Composition density credit mismatch")
    require(exact(composition["strict_upper_bound"]) == threshold, "Composition strict bound mismatch")
    require(objective < threshold < Q(63, 100), "Final exact inequality fails")
    require(Q((objective*10**12).numerator//(objective*10**12).denominator+1, 10**12) == threshold,
            "Strict bound is not the expected upward rounding")
    output.mkdir(parents=True, exist_ok=True)
    radial = {"status": "SEVEN_CLOSED_RADIAL_SIGNS_INDEPENDENTLY_CERTIFIED",
              "proof_method": "Rational Bernstein subdivision via de Casteljau identity",
              "source_sha256": hashes, "regions": signs,
              "nodes": sum(p["nodes"] for p in signs),
              "terminal_boxes": sum(p["terminal_boxes"] for p in signs)}
    vectors = {"status": "EXACT_RECONSTRUCTED_ORIGIN_DATA", "source_sha256": hashes,
               "matrix_definition": "M=(sum_j w_j q_j q_j^T)/common_denominator; q_j=ray_j.vector[1:] in kernel.blocks[0]",
               "dimension": n, "common_denominator": str(denominator),
               "degree_zero_positive_rays": len(block["rays"]),
               "u": list(map(str, u)), "v": list(map(str, v)), "kappa": str(kappa),
               "a": list(map(str, a)), "b": list(map(str, b)),
               "u_minus_b": list(map(str, unew)), "v_minus_a": list(map(str, vnew)),
               "new_psd_kappa": str(kappa-gain), "g": str(gain),
               "physical_pair_origin_v_unchanged": origin["pair_origin_v"],
               "physical_origin_warning": "pair_origin_v, v, and functional_origin_v have different roles. The PSD shift v->v-a is restored by the added pure-linear energy +2 sum w; no physical origin field is replaced.",
               "physical_pair_convention_from_source": origin["pair_convention"]}
    objective_data = {"status": "EXACT_FINAL_OBJECTIVE_ARITHMETIC_VERIFIED",
                      "scope": "Graph capacities and the global sign proof are separate obligations.",
                      "base_objective": str(base), "base_objective_decimal": dec(base),
                      "graph_capacity_order": names, "capacities": capacities,
                      "weights": list(map(str, weights)),
                      "graph_density_charge": str(graph_charge),
                      "raw_objective": str(raw), "raw_objective_decimal": dec(raw),
                      "origin_credit_g": str(gain), "density_credit": str(gain/8),
                      "density_credit_decimal": dec(gain/8),
                      "final_objective": str(objective), "final_objective_decimal": dec(objective),
                      "strict_upper_bound": str(threshold), "strict_upper_bound_decimal": dec(threshold),
                      "threshold_minus_final_objective": str(threshold-objective),
                      "threshold_minus_final_objective_decimal": dec(threshold-objective),
                      "final_objective_strictly_below_threshold": True,
                      "threshold_strictly_below_63_100": True,
                      "source_sha256": hashes}
    record = {"status": "PASS_INDEPENDENT_EXACT_ORIGIN_CREDIT_AND_OBJECTIVE_AUDIT",
              "scope": "Group 3 subcertificate only; no global branch verification performed by this checker.",
              "immutable_source_files": hashes,
              "matrix_dimension": n, "positive_gram_rays": len(block["rays"]),
              "range_equations_Mu_v": n, "range_equations_Mb_a": n,
              "shifted_range_equations": n, "radial_closed_regions": len(signs),
              "radial_nodes": radial["nodes"], "radial_terminal_boxes": radial["terminal_boxes"],
              "physical_origin_coupling_preserved_by_exact_cancellation": True,
              "pair_and_cyclic_expressions_preserved": True,
              "final_objective_decimal": dec(objective),
              "strict_upper_bound": str(threshold),
              "all_math_is_exact_rational": True,
              "uses_stored_radial_success_flags": False,
              "imports_original_verification_code": False,
              "dependencies": ["Python 3.9+ standard library"]}
    for name, result in [("ORIGIN_CREDIT_AUDIT.json", record),
                         ("origin_credit_vectors.json", vectors),
                         ("radial_bernstein_proof.json", radial),
                         ("final_objective.json", objective_data)]:
        (output/name).write_text(json.dumps(result, indent=2)+"\n")
    require(hashes == {p: sha(archive/p) for p in paths.values()}, "Original input files changed during audit")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    audit(args.archive, args.output)
