"""Check dated literature comparisons by exact rational inequalities.

PRIORITY.md and remaining_dimensions_review.json record the source review.
Vardy's decimal prefixes require intervals: prove the gain even at the upper
end of each reported prefix. Numerical checks do not automate bibliography.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import json
import argparse

if not __debug__:
    raise RuntimeError('Run without -O.')
P = Path(__file__).resolve().parent
SCALE = 10**6


def dec(x):
    x = F(x)
    return Decimal(x.numerator)/Decimal(x.denominator)


def audit():
    review = json.loads((P/'literature_review.json').read_text())
    data = json.loads((P.parent/'range48_96_20261115/results.json').read_text())
    by_dimension = {a['dimension']: a for a in data['rows']}
    prefixes = review['historical_guarantee_log2_prefixes']
    rows = []
    with localcontext() as ctx:
        ctx.prec = 70
        for n in review['registered_dimensions']:
            a = by_dimension[n]
            ours = F(a['best_checked']['density_squared'])
            prior = F(a['comparison']['density_squared'])
            assert ours > prior
            result = {'dimension': n, 'our_density_squared': str(ours),
                      'our_method': a['best_checked'],
                      'named_external_comparison': a['comparison'],
                      'comparison_verified_exactly': True,
                      'literature_status': 'best_known_lower_bound_in_dated_review',
                      'review_date': review['review_date']}
            if 'density_squared_interval' in a['comparison']:
                lower,upper = map(F,a['comparison']['density_squared_interval'])
                assert prior == lower and ours > upper >= lower
                ratio_lo,ratio_hi = ours/upper,ours/lower
                gain_lo = 100*(dec(ratio_lo).sqrt()-1)
                gain_hi = 100*(dec(ratio_hi).sqrt()-1)
                lo,hi = F(int(gain_lo*SCALE),SCALE),F(int(gain_hi*SCALE)+1,SCALE)
                assert (1+lo/100)**2 <= ratio_lo <= ratio_hi < (1+hi/100)**2
                result.update(comparison_kind='derived_classical_guarantee_interval',
                              comparison_density_squared_interval=[str(lower),str(upper)],
                              gain_percent_interval=[str(dec(lo)),str(dec(hi))],
                              gain_percent_display=f'>{dec(lo):.6f}%',
                              gain_uses_upper_endpoint_of_full_external_theorem_evaluation=True)
            elif str(n) in prefixes:
                lower = F(prefixes[str(n)])
                upper = lower+F(1, 100)
                gain_lo = 100*(dec(ours).sqrt()/Decimal(2)**dec(upper)-1)
                gain_hi = 100*(dec(ours).sqrt()/Decimal(2)**dec(lower)-1)
                k_lo = int(gain_lo*SCALE)
                k_hi = int(gain_hi*SCALE)+1
                lo, hi = F(k_lo, SCALE), F(k_hi, SCALE)
                assert lo > 0
                # Bounds on the reported guarantee B, not on the unknown
                # actual density: 2^lower <= B < 2^upper.
                e_lo, e_hi = 2*lower, 2*upper
                assert ours**e_hi.denominator > 2**e_hi.numerator*(1+lo/100)**(2*e_hi.denominator)
                assert ours**e_lo.denominator < 2**e_lo.numerator*(1+hi/100)**(2*e_lo.denominator)
                result.update(comparison_kind='published_decimal_prefix',
                              reported_guarantee_log2_interval=[str(lower), str(upper)],
                              reported_guarantee_log2_interval_lower_inclusive=True,
                              reported_guarantee_log2_interval_upper_inclusive=False,
                              gain_percent_interval=[str(dec(lo)), str(dec(hi))],
                              gain_percent_interval_endpoints_strict=True,
                              gain_percent_display=f'>{dec(lo):.6f}%',
                              gain_percent_range_display=f'{dec(lo):.6f}–{dec(hi):.6f}%',
                              exact_integer_power_comparisons_verified=True)
            else:
                ratio2 = ours/prior
                gain = 100*(dec(ratio2).sqrt()-1)
                k = int(gain*SCALE)
                lo, hi = F(k, SCALE), F(k+1, SCALE)
                assert (1+lo/100)**2 <= ratio2 < (1+hi/100)**2
                result.update(comparison_kind='exact_squared_density',
                              comparison_density_squared=str(prior),
                              ratio_squared_exact=str(ratio2), gain_percent=str(gain),
                              gain_percent_interval=[str(dec(lo)), str(dec(hi))],
                              gain_percent_display=f'{gain:.6f}%')
            rows.append(result)
        # The printed Prop3.2.22 construction has a norm6 vector. The
        # period-repaired lamination is a valid derived comparison only.
        d54 = F(5**54, 2**54*3**25)
        d55 = d54*F(9, 4)
        reported = d54*F(15,4)
        tilted = json.loads((P.parent/'tilted55_20261116/results.json').read_text())
        assert F(by_dimension[55]['comparison']['density_squared']) == d55
        assert F(by_dimension[55]['best_checked']['density_squared']) == F(tilted['density_squared'])>d55
        assert F(tilted['reported_external_claim_density_squared'])==reported
        assert not tilted['world_record_claimed']
        old = json.loads((P/'sections_results.json').read_text())['improved_controls']
        old55 = F(next(a for a in old if a['dimension'] == 55)['density_squared'])
        assert d55 > old55 and 55 not in review['registered_dimensions']
        d55_decimal = str(dec(d55).sqrt())
    remaining = json.loads((P/'remaining_dimensions_review.json').read_text())
    assert remaining['all_best_known_construction_dimensions'] == review['registered_dimensions']
    groups = [set(review['registered_dimensions']), set(remaining['classical_derived_controls']),
              set(remaining['dimensions_without_registered_record_claims'])]
    assert all(not a & b for i, a in enumerate(groups) for b in groups[i+1:])
    assert set.union(*groups) == set(range(48, 97))
    return {'review': review, 'rows': rows, 'all_comparisons_verified_exactly': True,
            'all_dimensions_classified': True,
            'Leterrier55_source_correction': {
                'source': remaining['sources']['Leterrier2024'],
                'density_squared': str(d55), 'center_density': d55_decimal,
                'comparison_kind':'derived_period_repair_not_the_printed_density',
                'printed_density_squared_not_used_as_guarantee':str(reported),
                'explicit_counterexample_norm_squared':6,
                'new_certified_55D_density_squared':tilted['density_squared'],
                'world_record_claimed':False,
                'earlier_repository_control_squared': str(old55),
                'repaired_lamination_exceeds_earlier_repository_control': True}}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--write-result', action='store_true')
    args = ap.parse_args()
    result = audit()
    path = P/'literature_comparisons.json'
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert json.loads(path.read_text()) == result, 'Stored literature comparison is stale.'
    for a in result['rows']:
        print(a['dimension'], a['gain_percent_display'], 'comparison verified')
    print('55D: repaired lamination comparison', result['Leterrier55_source_correction']['center_density'],
          '; printed claim disputed; no new world-record claim')
