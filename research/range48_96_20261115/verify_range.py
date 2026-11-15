"""Exact lower-bound survey of EVERY dimension 48--96; not a SOTA catalogue.

New computations: P48 compatible-half-coset counting and a broad colored
torus scan. Published/derived comparisons and repository results are separate.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
import json
import sys

if not __debug__:
    raise RuntimeError("Run without -O.")
P = Path(__file__).resolve().parent
ROOT = P.parents[1]
sys.path.insert(0, str(P.parent/'cross_dimension_lower_20261115'))
from verify_counting import power, mul, sigma, ORDER, determinant, root_gram, decsqrt
from verify_to96 import pi_interval, hermite_upper, section_floor, exceptional
sys.path.insert(0, str(P.parent/'priority_resolution_20261115'))
from verify_torus_extensions import color_data


def read(path):
    return json.loads((ROOT/path).read_text())


def row(n, name, d2, evidence, kind):
    return {'dimension': n, 'name': name, 'density_squared': str(F(d2)),
            'evidence': evidence, 'kind': kind}


def parent_inputs():
    mw = 'https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf'
    bcs = 'https://neilsloane.com/doc/Me103.pdf'
    return [row(48, 'Leech--Sloane P48q', F(3, 2)**48, bcs, 'published_lattice'),
            row(54, 'Leterrier54, Corollary 1.4 (2022)', F(5**54, 2**54*3**25), 'https://doi.org/10.1007/s40993-022-00321-0', 'published_lattice_guarantee'),
            row(56, 'Shioda MW56', F(5**56, 2**56*29**8), mw, 'published_lattice'),
            row(64, 'Nebe Ne64', F(3**32), 'https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html', 'published_lattice'),
            row(72, 'Nebe Gamma72', F(2**72), 'https://arxiv.org/abs/1008.2862', 'published_lattice'),
            row(80, 'Shioda MW80', F(7**80, 2**80*41**12), mw, 'published_lattice'),
            row(88, 'Classical BCS88', F(2**88), bcs, 'derived_classical_lattice'),
            row(92, 'Shioda MW92', F(2**184, 47**14), mw, 'published_lattice'),
            row(96, 'BCS96, Table 1 (1982)', F(3**96, 2**48), bcs, 'published_lattice')]


def p48_counting():
    e4 = [1]+[240*sigma(j, 3) for j in range(1, ORDER+1)]
    e6 = [1]+[-504*sigma(j, 5) for j in range(1, ORDER+1)]
    a, b = power(e4, 3), power(e6, 2)
    delta = [(x-y)//1728 for x, y in zip(a, b)]
    assert all(x-y == 1728*z for x, y, z in zip(a, b, delta))
    # Independent Delta expansion via the eta product.
    eta = [1]+[0]*ORDER
    from math import comb
    for j in range(1, ORDER+1):
        f = [0]*(ORDER+1)
        for k in range(min(24, ORDER//j)+1):
            f[j*k] = (-1)**k*comb(24, k)
        eta = mul(eta, f)
    assert delta == [0]+eta[:ORDER]
    terms = [mul(power(e4, 6-3*j), power(delta, j)) for j in range(3)]
    t, coefficients = terms[0][:], [1]
    for j in (1, 2):
        assert terms[j][j] == 1
        c = -t[j]
        coefficients.append(c)
        t = [x+c*y for x, y in zip(t, terms[j])]
    assert coefficients == [1, -1440, 125280]
    assert t[:4] == [1, 0, 0, 52416000]
    assert all(x >= 0 and x % 2 == 0 for x in t[3:])
    B0, B1 = t[4]//2, (t[3]+t[5])//2
    population = 2**47-2**23-1
    out, attempts = [], []
    for r in range(1, 9):
        Q = root_gram(r) if r <= 5 else exceptional(r)
        assert all(determinant([v[:j] for v in Q[:j]]) > 0 for j in range(1, r+1))
        counts = [0, 0]
        # Explicit embedding of q_R in r hyperbolic planes. Its e-coordinates
        # are the coefficient vector, which proves injectivity.
        w = []
        for i in range(r):
            e, f = 1 << i, 1 << i
            for j in range(i):
                f |= (Q[i][j] % 2) << j
            w.append((e, f))
        seen = set()
        for z in product((0, 1), repeat=r):
            e = f = 0
            for bit, (a, b) in zip(z, w):
                if bit:
                    e ^= a
                    f ^= b
            assert e not in seen
            seen.add(e)
            q = sum(z[i]*Q[i][j]*z[j] for i in range(r) for j in range(r))//2 % 2
            assert (e & f).bit_count() % 2 == q
            if any(z):
                counts[q] += 1
        expectation = F(counts[0]*B0+counts[1]*B1, population)
        attempts.append({'rank': r, 'type_counts': counts, 'expectation_upper': str(expectation), 'passes': expectation < 1})
        if expectation >= 1:
            continue
        G = [[F(5, 4)*x for x in v] for v in Q]
        # Nonsingular class: horizontal >=7/2, vertical >=5/2.
        # Nonzero singular class: horizontal >=3, vertical >=5.
        assert F(7, 2)+F(5, 2) == 6 and 3+5 >= 6 and 4*F(5, 2) >= 6
        d2 = F(3, 2)**48*F(6, 5)**r/determinant(Q)
        assert d2 == F(3, 2)**(48+r)/determinant(G)
        a = row(48+r, f'P48 compatible half-cosets, rank {r}', d2,
                'research/range48_96_20261115/P48_PROOF.md', 'new_existence_proof')
        a.update(root_gram=Q, vertical_gram=[[str(x) for x in v] for v in G],
                 expectation_upper=str(expectation), minimum_squared=6, explicit_shifts=False)
        out.append(a)
    assert [a['dimension'] for a in out] == list(range(49, 55))
    explicit = read('research/range48_96_20261115/explicit_results.json')
    assert explicit['status'] == 'verified_explicit_49_to_52_extensions'
    assert explicit['positive_cosets_completed'] == 15
    for a, b in zip(out[:4], explicit['rows']):
        assert a['dimension'] == b['dimension']
        assert a['density_squared'] == b['center_density_squared_exact']
        assert b['explicit_horizontal_shifts'] and b['minimum_squared'] == 6
        a['explicit_shifts'] = True
        a['kind'] = 'new_counting_proof_and_explicit_lattice'
    return {'theta_coefficients': coefficients, 'shell_counts': {str(2*j): t[j] for j in range(3, 7)},
            'bad_singular_upper': B0, 'bad_nonsingular_upper': B1,
            'population_lower': population, 'attempts': attempts, 'rows': out}


def torus_scan(parents):
    _, hi = pi_interval()
    out = []
    for parent in parents:
        n, d2 = parent['dimension'], F(parent['density_squared'])
        assert n % 2 == 0
        Vhi = hi**(n//2)/factorial(n//2)
        for r in range(1, min(12, 96-n)+1):
            if r <= 8:
                Q = root_gram(r) if r <= 5 else exceptional(r)
            else:
                A, B = exceptional(8), root_gram(r-8)
                Q = [v+[0]*len(B) for v in A]+[[0]*8+v for v in B]
            assert all(determinant([v[:j] for v in Q[:j]]) > 0 for j in range(1, r+1))
            schemes = [('mod2', 2**r, F(4), Q)]
            if r <= 4:
                Qc, C, N, lam = color_data(r)
                schemes.append(('small_period', N, lam, Qc))
            # A_r / (root-color kernel): all norm-2 roots change color.
            Qa = root_gram(r, 'A')
            schemes.append(('A_root_coloring', r+1, F(4 if r == 1 else 3 if r == 2 else 2), Qa))
            for scheme, N, lam, Qs in schemes:
                def fits(d):
                    return (N-1)**2*d2*Vhi**2*d**n < 1
                left, right = 0, 4_000_000
                assert not fits(F(right, 10**6))
                while left+1 < right:
                    mid = (left+right)//2
                    if fits(F(mid, 10**6)):
                        left = mid
                    else:
                        right = mid
                d = F(left, 10**6)
                assert fits(d) and not fits(d+F(1, 10**6))
                t = max(4-d, 4/lam)
                assert d+t >= 4 and lam*t >= 4
                result = d2*F(2, 1)**r/t**r/determinant(Qs)
                a = row(n+r, f'{parent["name"]} torus + {r}D {scheme}', result,
                        'research/range48_96_20261115/README.md', 'new_torus_calculation')
                a.update(parent_dimension=n, parent_density_squared=str(d2), colors=N,
                         vertical_root_gram=Qs, same_color_ratio=str(lam),
                         torus_squared_separation=str(d), vertical_squared_minimum=str(t),
                         strict_volume_inequality_checked=True, explicit_shifts=False)
                out.append(a)
    return out


def baseline(parents):
    lo, _ = pi_interval()
    out = list(parents)
    for p in parents:
        d2 = F(p['density_squared'])
        for n in range(p['dimension'], 48, -1):
            c = hermite_upper(n, lo)
            d2 = section_floor(d2, n, c)
            if d2 > 0:
                out.append(row(n-1, f'{p["name"]} / Blichfeldt section to {n-1}D', d2,
                               p['evidence'], 'derived_classical_section'))
            else:
                break
    classical = read('research/priority_resolution_20261115/classical_bcs88_sections.json')
    for a in classical['rows']:
        normal = {84:'D4',85:'A3',86:'A2',87:'A1'}[a['dimension']]
        out.append(row(a['dimension'], f'BCS88 section perpendicular to sqrt(2) {normal}', a['classical_center_density_squared_exact'],
                       'research/priority_resolution_20261115/CLASSICAL_SECTIONS.md', 'derived_classical_section'))
    for a in read('research/cross_dimension_lower_20261115/l80_results.json')['rows']:
        n = a['dimension']
        d2 = F(2**n)/F(a['normal_gram_determinant_exact'])
        out.append(row(n, f'Bachoc--Nebe L80 bare section, normal rank {80-n}', d2,
                       'research/cross_dimension_lower_20261115/l80_results.json', 'derived_parent_section_with_repository_embedding'))
    # Every norm-8 vector of the published unimodular L80 is primitive.
    # Its perpendicular gives the classical codimension-one section exactly.
    out.append(row(79, 'Bachoc--Nebe L80 perpendicular to a minimum vector', F(2**76),
                   'https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf', 'derived_classical_section'))
    sections = read('research/novel_routes_20261113/audits/gamma72_leech_sections_verified.json')
    for a in sections['sections']:
        out.append(row(a['dimension'], 'Classical Gamma72/Leech-block section', a['center_density_squared'],
                       'research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py', 'derived_classical_section'))
    out.append(row(59, 'Classical Gamma72/kappa13 section', F(2**44, 3**5),
                   'constructions/dim59/verify_literature_section.py', 'derived_classical_section'))
    # Printed prefixes are LOWER guarantees, not upper bounds on Elkies's
    # actual lattice densities. Exact integer powers certify each floor.
    for n, prefix in {57:'15.37',58:'16.46',59:'17.75',60:'19.04',61:'20.54',62:'22.13',63:'23.92'}.items():
        e = 2*F(prefix)
        target = 2**e.numerator
        low, high = 0, 1 << ((e.numerator+e.denominator-1)//e.denominator)
        while low+1 < high:
            mid = (low+high)//2
            if mid**e.denominator <= target:
                low = mid
            else:
                high = mid
        assert low**e.denominator <= target < (low+1)**e.denominator
        out.append(row(n, f'Elkies E{n}, Vardy1999 printed log2 prefix {prefix}', low,
                       'https://doi.org/10.1090/S0002-9947-99-02169-8', 'reported_historical_guarantee'))
    return close_products(out)


SMALL = {1:F(1,4),2:F(1,12),3:F(1,32),4:F(1,64),5:F(1,128),6:F(1,192),7:F(1,256),8:F(1,256),
         9:F(1,512),10:F(1,768),11:F(1,972),12:F(1,729),13:F(1,972),14:F(1,768),15:F(1,512),
         16:F(1,256),17:F(1,256),18:F(1,192),19:F(1,128),20:F(1,64),21:F(1,32),22:F(1,12),23:F(1,4),24:F(1)}


def close_products(rows):
    best = {}
    for a in rows:
        n = a['dimension']
        if 48 <= n <= 96 and (n not in best or F(a['density_squared']) > F(best[n]['density_squared'])):
            best[n] = a
    for n in range(48, 97):
        if n not in best:
            continue
        a = best[n]
        for k, s in SMALL.items():
            if n+k > 96:
                continue
            d2 = F(a['density_squared'])*s
            if n+k not in best or d2 > F(best[n+k]['density_squared']):
                best[n+k] = row(n+k, a['name']+f' x classical {k}D lattice', d2,
                                 a['evidence'], 'derived_product')
    return best


def prior_repo():
    rows = []
    for path, key in [('research/cross_dimension_lower_20261115/counting_results.json','rows'),
                      ('research/cross_dimension_lower_20261115/counting80_results.json','rows'),
                      ('research/cross_dimension_lower_20261115/l80_results.json','rows'),
                      ('research/priority_resolution_20261115/torus93_96_results.json','rows'),
                      ('research/novel_routes_20261113/higher_parent/gamma72_verified.json','certificates')]:
        for a in read(path)[key]:
            n = a['dimension']
            if 'counting80' in path:
                name = f'L80 compatible half-cosets, rank {n-80}'
            elif 'counting_results' in path:
                name = f'Gamma72 compatible half-cosets, rank {n-72}'
            elif 'l80_results' in path:
                name = f'L80 antipode, {a["cosets"]} translates'
            elif 'torus93' in path:
                name = f'MW92 periodic layers, {a["colors"]} colors'
            else:
                name = f'Gamma72 antipode, {a["cosets"]} translates'
            rows.append(row(n, name, a['center_density_squared_exact'], path, 'prior_repository_result'))
    for a in read('research/novel_routes_20261113/dual64/constructed.json'):
        rows.append(row(a['dimension'], f'Ne64 antipode, {a["cosets"]} translates', a['delta_squared'], 'research/novel_routes_20261113/dual64/constructed.json', 'prior_repository_result'))
    for a in read('research/novel_routes_20261113/audits/propagated_benchmarks.json')['all_finite_family_bounds']:
        rows.append(row(a['dimension'], a['name'], a['center_density_squared'], 'research/novel_routes_20261113/audits/propagated_benchmarks.json', 'prior_repository_derived_control'))
    for path, key in [('constructions/dim59/result.json','center_density_squared_exact'),
                      ('research/novel_routes_20261113/higher_parent/ne65_verified.json','center_density_squared')]:
        a = read(path)
        name = 'Ne64 antipode, six translates' if a['dimension'] == 59 else 'Explicit Ne65, certified depth 26/3'
        rows.append(row(a['dimension'], name, a[key], path, 'prior_repository_result'))
    return rows


def audit():
    parents = parent_inputs()
    review=read('research/covering_and_deep_cosets_20261115/literature_review.json')
    external = baseline(parents)
    previous = close_products(list(external.values())+prior_repo())
    counting = p48_counting()
    torus = torus_scan(parents)
    additional=[]
    path='research/covering_and_deep_cosets_20261115/deeper_results.json'
    deeper=read(path)
    assert deeper['status']=='verified_deeper_explicit_49_to_52'
    for a in deeper['rows']:
        additional.append(row(a['dimension'],'Deeper P48 cosets with fitted vertical metric',
                              a['density_squared'],path,'new_explicit_lattice'))
    path='research/covering_and_deep_cosets_20261115/covering_results.json'
    covering=read(path)
    assert covering['all_arithmetic_inequalities_exact']
    for a in covering['improved_dimensions']:
        additional.append(row(a['dimension'],a['name'],a['density_squared'],path,'new_covering_existence_proof'))
    path='research/covering_and_deep_cosets_20261115/sections_results.json'
    sections=read(path)
    assert sections['all_density_and_root_inequalities_exact']
    for a in sections['improved_controls']:
        additional.append(row(a['dimension'],a['name'],a['density_squared'],path,'new_derived_section_control'))
    new = close_products(list(previous.values())+counting['rows']+torus+additional)
    assert set(new) == set(previous) == set(external) == set(range(48, 97))
    results = []
    for n in range(48, 97):
        a, old, ext = new[n], previous[n], external[n]
        d2, prev2, ext2 = map(lambda x:F(x['density_squared']), (a, old, ext))
        results.append({'dimension':n, 'best_checked':a, 'center_density':decsqrt(d2),
                        'comparison':ext, 'comparison_center_density':decsqrt(ext2),
                        'previous_checked':old, 'previous_center_density':decsqrt(prev2),
                        'ratio_to_comparison':decsqrt(d2/ext2),
                        'ratio_to_previous':decsqrt(d2/prev2),
                        'improves_previous_checked_envelope':d2>prev2,
                        'literature_status':('best_known_lower_bound_in_dated_review'
                            if n in review['registered_dimensions'] and d2>ext2
                            else 'improves_named_comparison' if d2>ext2 else 'literature_input_or_derived_control'),
                        'literature_review_date':review['review_date']})
    return {'dimensions':list(range(48,97)), 'maximum_parent_dimension':96,
            'comparison_scope':'Published lower-bound guarantees and named derived classical controls; best-known claims use the dated source review in covering_and_deep_cosets_20261115/PRIORITY.md',
            'parent_inputs':parents, 'p48_counting':counting, 'torus_candidates':torus,
            'subsequent_audited_inputs':additional,
            'rows':results, 'all_inequalities_exact':True}


def render_table(result):
    def link(a):
        url = a['evidence'] if '://' in a['evidence'] else '../../'+a['evidence']
        return f'[{a["name"]}]({url})'
    lines = ['# Every dimension from 48D through 96D', '',
             'Generated by `verify_range.py`. All densities are approximate **center densities** of proved lower bounds, not occupied fractions.', '',
             'Best-known claims use the dated literature review in [PRIORITY.md](../covering_and_deep_cosets_20261115/PRIORITY.md). The comparison column names the strongest bound in the original published-input/Blichfeldt-section/product comparison family, before the covering and sharper-section updates. Derived numbers are explicitly distinguished from numbers printed by the source. We compare proved lower bounds; an uncomputed actual lattice minimum is not a prerequisite.', '',
             '| Dimension | Best bound in this audit | Construction / evidence | Specified comparison | Comparison density | Ratio | Literature status |',
             '| --- | ---: | --- | --- | ---: | ---: | --- |']
    for a in result['rows']:
        new = ('Best-known; review '+a['literature_review_date'] if a['literature_status']=='best_known_lower_bound_in_dated_review'
               else 'Improves named bound' if a['literature_status']=='improves_named_comparison' else 'Literature / derived control')
        comparison = a['comparison']
        kind = 'reported guarantee' if comparison['kind'] == 'reported_historical_guarantee' else 'published parent/guarantee' if comparison['kind'].startswith('published_') else 'derived'
        lines.append(f'| {a["dimension"]} | {Decimal(a["center_density"]):.9E} | {link(a["best_checked"])} | {link(comparison)} ({kind}) | {Decimal(a["comparison_center_density"]):.9E} | {Decimal(a["ratio_to_comparison"]):.6f} | {new} |')
    return '\n'.join(lines)+'\n'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--write-result', action='store_true')
    args = ap.parse_args()
    result = audit()
    p = P/'results.json'
    if args.write_result:
        p.write_text(json.dumps(result,indent=2)+'\n')
        (P/'TABLE.md').write_text(render_table(result))
    else:
        assert json.loads(p.read_text()) == result, 'Stored range report is stale.'
        assert (P/'TABLE.md').read_text() == render_table(result), 'Stored table is stale.'
    for a in result['rows']:
        print(a['dimension'],a['center_density'][:19],a['best_checked']['name'],
              'NEW' if a['improves_previous_checked_envelope'] else '')
