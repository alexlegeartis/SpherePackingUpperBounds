"""Exact source-to-period audit of the explicit 86D--91D MW92 antipodes.

Published parent minimum16 is inherited. Every curve equation, height,
Gram identity, full-basis identity, kernel and projected shift is checked.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import prod, factorial
import argparse, json, sys
import numpy as np
from flint import fmpz_mat, fmpq_mat, fq_default_ctx, fq_default_poly_ctx, fmpz_mod_poly_ctx

if not __debug__:
    raise RuntimeError('Run without -O.')
P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent/'cross_dimension_lower_20261115'))
from verify_to96 import pi_interval, section_floor
from verify_counting import decsqrt


def parent():
    source = json.loads((P/'unitary_sections.json').read_text())
    p = source['prime']
    assert p == 47 and source['field_nonsquare'] == -1
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(p)([1,0,1]), fq_type='FQ_NMOD')
    polyring = fq_default_poly_ctx(field)
    def element(z):
        assert isinstance(z, int) and 0 <= z < p*p
        return field([z % p, z//p])
    polys = []
    for row in source['unitaries']:
        aa,bb,cc,dd = map(element,row)
        assert aa*aa**p+bb*bb**p == 1
        assert cc*cc**p+dd*dd**p == 1
        assert aa*cc**p+bb*dd**p == 0
        x = -polyring([bb,aa])**16
        y = polyring([dd,cc])**24
        assert y*y == x*x*x+polyring([1]+[0]*47+[1])
        polys.append((x,y))
    count = len(polys)
    assert count == 180
    heights = [[0]*count for _ in range(count)]
    for i,(x,y) in enumerate(polys):
        heights[i][i] = 16
        for j,(u,v) in enumerate(polys[:i]):
            dx,dy = x-u,y-v
            if dx.is_zero() and dy.is_zero():
                height = 16
            elif dx.is_zero() and (y+v).is_zero():
                height = -16
            else:
                finite = dx.gcd(dy).degree()
                # At infinity use the minimal smooth chart x'=s^16*x,
                # y'=s^24*y. Identically zero differences impose no bound.
                orders = []
                if not dx.is_zero(): orders.append(16-dx.degree())
                if not dy.is_zero(): orders.append(24-dy.degree())
                infinite = min(orders)
                assert finite >= 0 and infinite >= 0
                height = 8-finite-infinite
            heights[i][j] = heights[j][i] = height
    with np.load(P/'parent_certificate.npz',allow_pickle=False) as data:
        assert set(data.files) == {'G','section_coordinates'}
        assert all(np.issubdtype(data[k].dtype,np.integer) for k in data.files)
        G,Y = [fmpz_mat(data[k].tolist()) for k in ('G','section_coordinates')]
    assert G.nrows() == G.ncols() == 92 and Y.nrows() == 180 and Y.ncols() == 92
    assert G == G.transpose() and all(G[i,i] % 2 == 0 for i in range(92))
    assert all(fmpz_mat([row[:k] for row in G.tolist()[:k]]).det() > 0 for k in range(1,93))
    assert Y*G*Y.transpose() == fmpz_mat(heights)
    assert Y.hnf().tolist()[:92] == [[int(i == j) for j in range(92)] for i in range(92)]
    assert G.det() == 47**14
    return G, {'explicit_polynomial_sections':180,'curve_equations_checked':180,
               'pairwise_heights_checked_by_direct_polynomial_gcd':180*179//2,
               'coordinates_generate_the_full_integer_module':True,
               'gram_determinant':str(G.det()),'minimum_squared_inherited':16,
               'source':'Schutt-Shioda Theorems6.23 and13.35; Section13.4.2',
               'no_new_parent_minimum_enumeration_claimed':True}


def rogers_comparison_intervals():
    # R_n denotes the numerical Rogers simplex upper bound, not the unknown
    # optimum packing density. Jensen bounds R_n FROM BELOW. This yields
    # an UPPER endpoint for the section guarantee derived from that theorem.
    _,pi_hi = pi_interval()
    density_upper = F(2**184,47**14)
    lower_rows = json.loads((P.parent/'covering_and_deep_cosets_20261115/sections_results.json').read_text())['improved_controls']
    lower = {a['dimension']:F(a['density_squared']) for a in lower_rows}
    out = {}
    for n in (92,91,90):
        volume_hi = (pi_hi**(n//2)/factorial(n//2) if n%2 == 0 else
                     F(2**(n//2+1),prod(range(1,n+1,2)))*pi_hi**(n//2))
        simplex_squared_lower = F((n+1)**2)*F(n+1,2*(n+3))**n
        gamma_power_lower = 4**n*simplex_squared_lower/volume_hi**2
        with localcontext() as ctx:
            ctx.prec = 80
            value = Decimal(gamma_power_lower.numerator)/gamma_power_lower.denominator
            proposed = int((value.ln()/n).exp()*10**12)
        c = F(proposed,10**12)
        while c**n > gamma_power_lower: c -= F(1,10**12)
        assert 0 < c**n <= gamma_power_lower
        density_upper = section_floor(density_upper,n,c)+1
        assert lower[n-1] < density_upper
        out[n-1] = {'density_squared_lower':str(lower[n-1]),
                    'density_squared_upper':str(density_upper),
                    'rogers_hermite_formula_lower':str(c),
                    'upper_endpoint_is_not_asserted_as_an_achieved_packing':True}
    return out


def gain_interval(ratio_lo,ratio_hi=None):
    if ratio_hi is None: ratio_hi = ratio_lo
    with localcontext() as ctx:
        ctx.prec = 70
        a = 100*((Decimal(ratio_lo.numerator)/ratio_lo.denominator).sqrt()-1)
        b = 100*((Decimal(ratio_hi.numerator)/ratio_hi.denominator).sqrt()-1)
        lo,hi = F(int(a*10**6),10**6),F(int(b*10**6)+1,10**6)
        assert (1+lo/100)**2 <= ratio_lo <= ratio_hi < (1+hi/100)**2
        return {'interval':[str(lo),str(hi)],'display':f'>{Decimal(lo.numerator)/lo.denominator:.6f}%' if ratio_lo != ratio_hi else f'{a:.6f}%',
                'approximate_percent':str(a)}


def audit(export=None):
    G,parent_report = parent()
    source = json.loads((P/'sections.json').read_text())
    checkpoint = {a['dimension']:a for a in json.loads((P/'checkpoint.json').read_text())['rows']}
    rogers = rogers_comparison_intervals()
    invG = G.inv()
    rows = []
    if export: export.mkdir(parents=True,exist_ok=False)
    for a in source['rows']:
        n = a['dimension']; k = 92-n
        B = fmpz_mat(a['normal_labels']); S = fmpz_mat(a['cluster_labels'])
        assert B.nrows() == k and B.ncols() == 92 and S.ncols() == k
        assert len({tuple(row) for row in S.tolist()}) == S.nrows()
        H,T = B.transpose().hnf(transform=True)
        assert T*B.transpose() == H and abs(T.det()) == 1
        assert H.tolist() == [[int(i == j) for j in range(k)] for i in range(92)]
        lifts = fmpz_mat(T.tolist()[:k]); N = fmpz_mat(T.tolist()[k:])
        assert lifts*B.transpose() == fmpz_mat([[int(i == j) for j in range(k)] for i in range(k)])
        assert N*B.transpose() == fmpz_mat(n,k)
        K = fmpq_mat(B)*invG*fmpq_mat(B.transpose())
        period = N*G*N.transpose()
        assert period.det() == G.det()*K.det()
        shifts = fmpq_mat(S*lifts)-fmpq_mat(S)*K.inv()*fmpq_mat(B)*invG
        assert shifts*fmpq_mat(B.transpose()) == fmpq_mat(S.nrows(),k)
        shift_coordinates = shifts*fmpq_mat(G*N.transpose())*period.inv()
        assert shift_coordinates*fmpq_mat(N) == shifts
        beta = F(0)
        for i in range(S.nrows()):
            for j in range(i):
                d = fmpq_mat([[int(S[i,t]-S[j,t]) for t in range(k)]])
                beta = max(beta,F(str((d*K.inv()*d.transpose())[0,0])))
        minimum = 16-beta
        assert 0 <= beta < 16
        d2 = F(S.nrows()**2)*(minimum/4)**n/F(int(period.det()))
        previous = F(checkpoint[n]['best_checked']['density_squared'])
        assert d2 > previous
        if n >= 89:
            comparison = rogers[n]
            ext_lo,ext_hi = map(F,(comparison['density_squared_lower'],comparison['density_squared_upper']))
            ext_name = 'MW92 section using the full Rogers simplex theorem (derived interval)'
        else:
            ext_lo = ext_hi = {86:F(2**84,3),87:F(2**85),88:F(2**88)}[n]
            comparison = {'density_squared_lower':str(ext_lo),'density_squared_upper':str(ext_hi)}
            ext_name = {86:'BCS88 section perpendicular to sqrt(2) A2',87:'BCS88 section perpendicular to sqrt(2) A1',88:'BCS88 E8/code construction'}[n]
        assert d2 > ext_hi
        row = {'dimension':n,'name':f'Explicit MW92 antipode, {S.nrows()} translates',
               'kind':'new_explicit_periodic_packing','density_squared':str(d2),'center_density':decsqrt(d2),
               'parent_dimension':92,'normal_rank':k,'normal_gram':[[str(x) for x in row] for row in K.tolist()],
               'normal_determinant':str(K.det()),'period_determinant':str(period.det()),
               'cosets':S.nrows(),'cluster_squared_diameter':str(beta),'minimum_squared_guarantee':str(minimum),
               'incremental_gain':gain_interval(d2/previous),
               'external_comparison_name':ext_name,'external_comparison':comparison,
               'external_gain':gain_interval(d2/ext_hi,d2/ext_lo),
               'complete_integer_kernel_checked':True,'all_projected_shifts_checked_exactly':True}
        rows.append(row)
        if export:
            (export/f'dim{n}.json').write_text(json.dumps({'period_basis_in_parent':N.tolist(),
                'period_gram':period.tolist(),'shifts_in_period_coordinates':[[str(x) for x in row] for row in shift_coordinates.tolist()],
                'shifts_in_parent':[[str(x) for x in row] for row in shifts.tolist()]},indent=2,default=int)+'\n')
    assert [a['dimension'] for a in rows] == list(range(86,92))
    return {'status':'six_new_explicit_dimension_improvements','maximum_parent_or_target_dimension':92,
            'parent_audit':parent_report,'rows':rows,'all_source_and_period_checks_exact':True,
            'checkpoint_commit':'a30fecabce649468f2d7f3ac1fc94c4b54e4c33e'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-result',action='store_true')
    parser.add_argument('--export',type=Path)
    args = parser.parse_args()
    result = audit(args.export)
    if args.write_result: (P/'results.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
    else: assert json.loads((P/'results.json').read_text()) == result
    for row in result['rows']:
        print(row['dimension'],row['center_density'],'previous',row['incremental_gain']['display'],'external',row['external_gain']['display'])
