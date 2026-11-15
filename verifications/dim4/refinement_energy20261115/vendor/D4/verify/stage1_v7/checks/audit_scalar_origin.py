"""Supplement fixed V7 replay with physical-row, objective and rejection checks.

Reads the supplied records only. No optimization, discovery, or global sign run.
"""
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
import argparse
import importlib.util
import json
import math
import sys

sys.dont_write_bytecode = True
# Stage 1 export: paths come from --record (a working copy of the record);
# the origin-credit module sits next to this file.
HERE = Path(__file__).resolve().parent
RECORD = None
SOURCE = None


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(value):
    require(type(value) in (int, str), 'Exact scalar requires integer or string')
    return Q(value)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dot(a, b):
    require(len(a) == len(b), 'Vector dimensions differ')
    return sum((x*y for x, y in zip(a, b)), Q(0))


def decimal(value):
    with localcontext() as context:
        context.prec = 70
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def laguerre_fourier(f):
    # R_k=k! L_k^(1); R_(k+1)=(2k+2-s)R_k-k(k+1)R_(k-1).
    rows = [[Q(1)], [Q(2), Q(-1)]]
    while len(rows) < len(f):
        k = len(rows)-1
        following = [Q(0)]*(k+2)
        for i, value in enumerate(rows[-1]):
            following[i] += (2*k+2)*value
            following[i+1] -= value
        for i, value in enumerate(rows[-2]):
            following[i] -= k*(k+1)*value
        rows.append(following)
    result = [Q(0)]*len(f)
    for multiplier, row in zip(f, rows):
        for i, value in enumerate(row):
            result[i] += multiplier*value
    return result


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def rational_remainder(a, b):
    result = list(map(Q, a))
    divisor = list(map(Q, b))
    require(divisor[-1] != 0, 'Zero divisor leading coefficient')
    while len(result) >= len(divisor) and result != [Q(0)]:
        shift = len(result)-len(divisor)
        factor = result[-1]/divisor[-1]
        for j, value in enumerate(divisor):
            result[j+shift] -= factor*value
        result = trim(result)
    return result


def prove_half_line(scalar, poly, left, wanted, check_sequence=False):
    sequence = scalar.sturm(scalar.integers(poly))
    require(len(sequence) >= 2 and sequence[-1] != [0], 'Degenerate Sturm sequence')
    identities = 0
    if check_sequence:
        derivative = [Q(i)*sequence[0][i] for i in range(1, len(sequence[0]))]
        ratio = derivative[-1]/sequence[1][-1]
        require(ratio > 0 and derivative == [ratio*x for x in sequence[1]],
                'Sturm derivative is not a positive rescaling')
        for a, b, following in zip(sequence, sequence[1:], sequence[2:]):
            negative = [-x for x in rational_remainder(a, b)]
            require(len(negative) == len(following), 'Sturm remainder degree differs')
            scale = negative[-1]/following[-1]
            require(scale > 0 and negative == [scale*x for x in following],
                    'Sturm recurrence is not a positive rescaling of negative remainder')
            identities += 1
        require(rational_remainder(sequence[-2], sequence[-1]) == [Q(0)],
                'Incomplete terminal Sturm remainder')
    a = [scalar.sign_at(p, left) for p in sequence]
    b = [scalar.sign_at(p, None) for p in sequence]
    va, vb = scalar.variations(a), scalar.variations(b)
    require(va == vb and a[0] == b[0] == wanted, 'Scalar half-line sign fails')
    return {'left': str(left), 'left_variations': va, 'infinity_variations': vb,
            'root_count': va-vb, 'polynomial_sign': wanted,
            'sturm_length': len(sequence), 'exact_rational_remainder_identities': identities,
            'terminal_remainder_zero': check_sequence}


def elevate_quartic(coefficients):
    require(len(coefficients) == 5, 'Quartic block dimension differs')
    return [sum((Q(math.comb(4, j)*math.comb(2, i-j), math.comb(6, i))*coefficients[j]
                 for j in range(max(0, i-2), min(4, i)+1)), Q(0)) for i in range(7)]


def check_functional(data, breaks):
    construction = data['primary_basis_construction']
    require(construction['kind'] == 'two_quartic_Bernstein_blocks_then_five_sorted_count_indicators',
            'Unexpected functional model')
    coefficients = list(map(exact, data['functional_origin_v']))
    require(len(coefficients) == 15, 'Functional origin dimension differs')
    radii = list(map(exact, construction['indicator_radii']))
    require(len(radii) == 5 and radii == breaks[1:6], 'Functional indicator radii differ')
    lifted = []
    for left, right in zip(breaks, breaks[1:]):
        row = elevate_quartic(coefficients[0:5] if right <= Q(3, 2) else coefficients[5:10])
        offset = sum((value for value, radius in zip(coefficients[10:], radii) if right <= radius), Q(0))
        lifted.extend(value+offset for value in row)
    require(lifted == list(map(exact, data['v'])), 'Functional origin differs from PSD physical basis')
    return lifted


def check_physical(data, breaks):
    n = 49
    shift = [Q(0)]*n
    counts = {}
    for field, coefficient_key in [('pointed_caps', 'M'), ('pointed_tents', 'capacity')]:
        counts[field] = len(data[field])
        for record in data[field]:
            weight = exact(record['weight'])
            capacity = exact(record[coefficient_key])
            w = list(map(exact, record['w']))
            require(len(w) == n and weight >= 0 and capacity >= 0, 'Malformed geometric linear term')
            for i, value in enumerate(w):
                shift[i] += capacity*weight*value/2
    counts['independent_antipodal_allocations'] = len(data['independent_antipodal_allocations'])
    for record in data['independent_antipodal_allocations']:
        require(exact(record['origin_constant']) == 0, 'Allocation has an uncharged origin')
        weight, w = exact(record['weight']), list(map(exact, record['w']))
        require(len(w) == n and weight >= 0, 'Malformed allocation linear term')
        for i, value in enumerate(w):
            shift[i] += weight*value/2
    counts['appended_spherical_local_energy_rays'] = len(data['appended_spherical_local_energy_rays'])
    for record in data['appended_spherical_local_energy_rays']:
        weight, constant = exact(record['weight']), exact(record['B'])
        require(weight >= 0 and exact(record['global_origin']) == 0,
                'Signed spherical component has uncharged origin')
        require(exact(record['global_origin_coupling_coefficient']) == constant/2,
                'Signed spherical coupling convention differs')
        radius = exact(record['R'])
        for region, (left, right) in enumerate(zip(breaks, breaks[1:])):
            require(right <= radius or left >= radius, 'Spherical support crosses physical piece')
            if right <= radius:
                for i in range(7*region, 7*(region+1)):
                    shift[i] += weight*constant/2
    counts['conditional_radial_rays'] = len(data['conditional_radial_rays'])
    for record in data['conditional_radial_rays']:
        require(record['kind'] == 'one_outer' and record['maximum'] == 24,
                'Unexpected conditional radial law')
        require(exact(record['pair_coefficient']) == 23,
                'Unexpected conditional pair coefficient')
        r, R = exact(record['inner_radius']), exact(record['outer_radius'])
        weight = exact(record['weight'])
        require(weight >= 0, 'Negative conditional radial weight')
        for region, (left, right) in enumerate(zip(breaks, breaks[1:])):
            require((right <= r or left >= r) and (right <= R or left >= R),
                    'Conditional radial support crosses physical piece')
            if left >= r and right <= R:
                for i in range(7*region, 7*(region+1)):
                    shift[i] += Q(23, 2)*weight
    v = list(map(exact, data['v']))
    pair = list(map(exact, data['pair_origin_v']))
    require(len(v) == len(pair) == n, 'Physical row dimensions differ')
    require([x+y for x, y in zip(v, shift)] == pair, 'Physical pair row does not equal PSD row plus all linear shifts')
    require(all(exact(ray['origin_value']) == 0 for ray in data['appended_broad19_nonorigin_psd_rays']),
            'Extra nonorigin PSD ray has a nonzero origin')
    forbidden = ['integer_quadratic', 'local_energy', 'count_gauge_products',
                 'product_count_gauges', 'tail_budget_caps']
    require(all(not data.get(name) for name in forbidden), 'Unaccounted local-energy representation')
    return shift, counts


def check_objective(data, scalar, graph, credit):
    f0 = exact(scalar['F'][0])
    kappa = exact(data['kappa'])
    nested = exact(data['nested_count_quadratic']['coefficients'][0])
    cuts = sum((exact(c['M'])*exact(c['weight']) for c in data.get('count_cuts', [])), Q(0))
    conditional = sum((exact(c['weight'])*exact(c['coefficients'][0])
                       for c in data['conditional_count_polynomials']), Q(0))
    base = (f0+kappa+nested+cuts+conditional)/8
    require(base == exact(data['objective']), 'Stored base objective does not match all origin constants')
    capacities = [862, 1102, 1358, 120, 604, 1880]
    weights = list(map(exact, graph['weights']))+[exact(graph[k]) for k in
              ('outer_graph_weight', 'narrow_graph_weight', 'moderate_graph_weight')]
    require(len(weights) == 6 and min(weights) >= 0, 'Graph weights invalid')
    charge = dot(capacities, weights)/8
    require(charge == Q(4769, 200000000), 'Graph charge differs from full V7 budget')
    gain = exact(credit['origin_credit'])
    require(exact(credit['density_credit']) == gain/8, 'Origin credit density normalization differs')
    final = base+charge-gain/8
    threshold = Q(78749932299, 125000000000)
    require(final < threshold, 'Exact final objective misses threshold')
    return {'F0': str(f0), 'PSD_kappa': str(kappa), 'nested_count_constant': str(nested),
            'count_cut_constant': str(cuts), 'conditional_count_constant': str(conditional),
            'base_local_origin_constant': str(kappa+nested+cuts+conditional),
            'all_graph_origin_constant': str(8*charge),
            'final_local_origin_constant': str(kappa+nested+cuts+conditional+8*charge-gain),
            'base_objective': str(base), 'graph_charge': str(charge),
            'density_credit': str(gain/8), 'final_objective': str(final),
            'final_objective_decimal': decimal(final), 'strict_threshold': str(threshold),
            'strict_gap': str(threshold-final), 'strict_gap_decimal': decimal(threshold-final)}


def main():
    global RECORD, SOURCE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RECORD = args.record.resolve()
    SOURCE = RECORD/'current/candidates/round9_six_graphs_degree18'
    require(not args.output.exists(), 'Output must be a new file')
    require(not args.output.resolve().is_relative_to(RECORD), 'Output lies inside supplied source')
    kernel, origin, scalar_data, graph = [json.loads((SOURCE/name).read_text()) for name in
                                        ('kernel.json', 'origin.json', 'scalar.json', 'graph_energies.json')]
    credit = json.loads((RECORD/'audit/round9_origin_credit/ORIGIN_CREDIT_VERIFIED.json').read_text())
    scalar = load('supplement_v7_scalar', SOURCE/'verify_scalar.py')
    radial = load('supplement_v7_radial', HERE/'check_origin_credit.py')
    f, g = [list(map(exact, scalar_data[name])) for name in ('F', 'G')]
    require(laguerre_fourier(f) == g and g[0] == 1, 'Independent Laguerre Fourier identity fails')
    scalar_proofs = {label: prove_half_line(scalar, poly, left, wanted, True)
                     for label, poly, left, wanted in
                     [('spatial_tail', f, Q(49, 8), -1), ('Fourier_half_line', g, Q(0), 1)]}
    breaks = list(map(exact, kernel['breaks']))
    check_functional(origin, breaks)
    shifts, component_counts = check_physical(origin, breaks)
    objective = check_objective(origin, scalar_data, graph, credit)
    controls = []

    def rejected(name, operation):
        try:
            operation()
        except (ValueError, AssertionError) as error:
            controls.append({'name': name, 'status': 'REJECTED',
                             'diagnostic': str(error)})
        else:
            raise ValueError('Negative control was accepted: '+name)

    bad = deepcopy(origin)
    bad['pair_origin_v'] = credit['new_psd_origin_v']
    rejected('substitute_shifted_PSD_row_for_physical_pair_row', lambda: check_physical(bad, breaks))
    missing = deepcopy(origin)
    missing['pointed_caps'] = missing['pointed_caps'][1:]
    rejected('omit_one_nonzero_cap_linear_shift', lambda: check_physical(missing, breaks))
    bad_objective = deepcopy(origin)
    bad_objective['objective'] = str(exact(bad_objective['objective'])-Q(1, 1000000))
    rejected('understate_stored_base_objective', lambda: check_objective(bad_objective, scalar_data, graph, credit))
    missing_graph = deepcopy(graph)
    missing_graph['moderate_graph_weight'] = '0'
    rejected('omit_moderate_graph_origin_charge', lambda: check_objective(origin, scalar_data, missing_graph, credit))
    bad_credit = deepcopy(credit)
    bad_credit['density_credit'] = str(2*exact(bad_credit['density_credit']))
    rejected('double_density_credit', lambda: check_objective(origin, scalar_data, graph, bad_credit))
    rejected('negative_radial_weight', lambda: radial.sign_certificate([Q(-1)]*7, Q(0), Q(1)))
    rejected('float_exact_input', lambda: radial.exact(0.5))
    rejected('Boolean_exact_input', lambda: radial.exact(True))
    rejected('Fourier_positive_endpoints_but_two_positive_roots',
             lambda: prove_half_line(scalar, [Q(1), Q(-3), Q(1)], Q(0), 1))
    rejected('Fourier_double_root_on_half_line',
             lambda: prove_half_line(scalar, [Q(1), Q(-2), Q(1)], Q(0), 1))
    rejected('tail_root_at_closed_cutoff',
             lambda: prove_half_line(scalar, [Q(49, 8), Q(-1)], Q(49, 8), -1))
    report = {'status': 'PASS_EXACT_PHYSICAL_ORIGIN_OBJECTIVE_STURM_AND_REJECTION_AUDIT',
              'scope': 'Fixed scalar/origin supplement only; global local-energy/cyclic prerequisites are separate.',
              'independent_fourier_recurrence': 'k! L_k^(1), exact three-term recurrence through k=49',
              'scalar_half_line_proofs': scalar_proofs,
              'functional15_to_Bernstein49_exact_equations': 49,
              'physical_pair_origin_exact_equations': 49,
              'linear_component_record_counts': component_counts,
              'physical_pair_origin_shift': list(map(str, shifts)),
              'objective_reconstruction': objective,
              'negative_controls': controls,
              'negative_control_count': len(controls)}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ('status', 'functional15_to_Bernstein49_exact_equations',
                    'physical_pair_origin_exact_equations', 'negative_control_count')}, indent=2))


if __name__ == '__main__':
    main()
