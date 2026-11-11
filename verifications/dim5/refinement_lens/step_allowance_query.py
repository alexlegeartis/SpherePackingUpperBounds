"""Exact conservative allowance queries for a previously verified record.

These helpers do not themselves certify its physical-budget inequalities.
Run construct_step_allowances.py --verify first and bind the physical payload.
Squared coordinates avoid introducing square-root rounding in far-box checks.
"""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
from fractions import Fraction as F


def base_squared_lower_bound(record,lo,hi):
    lo,hi=F(lo),F(hi)
    assert 0<=lo<=hi
    bins=record['bins'];boundaries=record['boundaries']
    first=F(boundaries[0]['base'])**2
    last=F(boundaries[-1]['base'])**2
    candidates=[]
    if lo<first or hi>last:
        candidates.append(F(0))
    for boundary in boundaries:
        p=F(boundary['base'])**2
        if lo<=p<=hi:
            candidates.append(F(boundary['eta']))
    for row in bins:
        a,b=F(row['lo'])**2,F(row['hi'])**2
        # Include every open bin meeting the closed query interval.
        if (lo==hi and a<lo<b) or (lo<hi and a<hi and b>lo):
            candidates.append(F(row['eta']))
    assert candidates,(lo,hi)
    return min(candidates)


def far_box_lower_bound(record,a2lo,a2hi,b2lo,b2hi,c2lo,c2hi):
    a2lo,a2hi,b2lo,b2hi=map(F,[a2lo,a2hi,b2lo,b2hi])
    assert 1<=a2lo<=a2hi and 1<=b2lo<=b2hi
    r2=F(record['short_radius'])**2
    if a2hi>r2 or b2hi>r2:
        # Split the radial box if it straddles the eligibility boundary and
        # retaining a positive allowance on the inside matters.
        return F(0)
    return base_squared_lower_bound(record,c2lo,c2hi)


def self_test(record):
    for row in record['bins']:
        a,b=F(row['lo'])**2,F(row['hi'])**2
        midpoint=(a+b)/2
        assert base_squared_lower_bound(record,midpoint,midpoint)==F(row['eta'])
    for boundary in record['boundaries']:
        x=F(boundary['base'])**2
        assert base_squared_lower_bound(record,x,x)==F(boundary['eta'])
    assert base_squared_lower_bound(record,0,100)==0
    assert far_box_lower_bound(record,1,2,1,1,F(19,10)**2,F(191,100)**2)==0
    return True


if __name__=='__main__':
    import json
    from pathlib import Path
    path=Path(__file__).with_name('step_allowances_005.json')
    data=json.loads(path.read_text())
    for record in data['dimensions']:
        assert self_test(record)
        print(record['dimension'],'exact query tests passed')
