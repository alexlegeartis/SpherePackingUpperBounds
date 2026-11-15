"""Check registered literature comparisons numerically; retain dated review inputs.

The bibliography review is recorded in literature_review.json and PRIORITY.md.
This checker certifies the density comparisons, not an automated web search.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json,argparse

if not __debug__:raise RuntimeError('Run without -O.')
P=Path(__file__).resolve().parent


def audit():
    review=json.loads((P/'literature_review.json').read_text())
    data=json.loads((P.parent/'range48_96_20261115/results.json').read_text())
    rows=[]
    for n in review['registered_dimensions']:
        a=next(a for a in data['rows']if a['dimension']==n)
        ours=F(a['best_checked']['density_squared']);prior=F(a['comparison']['density_squared'])
        assert ours>prior
        ratio2=ours/prior
        with localcontext()as ctx:
            ctx.prec=60
            ratio=(Decimal(ratio2.numerator)/Decimal(ratio2.denominator)).sqrt()
            gain=100*(ratio-1);k=int(gain*10**6)
        lo,hi=F(k,10**6),F(k+1,10**6)
        assert (1+lo/100)**2<=ratio2<(1+hi/100)**2
        rows.append({'dimension':n,'our_density_squared':str(ours),
                     'named_external_comparison':a['comparison'],
                     'comparison_density_squared':str(prior),'ratio_squared_exact':str(ratio2),
                     'gain_percent':str(gain),'gain_percent_interval':[str(lo),str(hi)],
                     'comparison_verified_exactly':True,
                     'literature_status':'best_known_lower_bound_in_dated_review',
                     'review_date':review['review_date']})
    return {'review':review,'rows':rows,'all_comparisons_verified_exactly':True}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args()
    result=audit();path=P/'literature_comparisons.json'
    if args.write_result:path.write_text(json.dumps(result,indent=2)+'\n')
    else:assert json.loads(path.read_text())==result,'Stored literature comparison is stale.'
    for a in result['rows']:print(a['dimension'],a['gain_percent'],'percent; comparison verified')
