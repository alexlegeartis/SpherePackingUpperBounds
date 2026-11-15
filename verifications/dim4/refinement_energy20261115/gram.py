"""Exact Schur credits from unused positive angular Gram blocks.

For M>=0 and Mw=v, beta=1/(v^T w) gives M-beta vv^T>=0.
The input M is independently rebuilt from the upstream positive integer rays.
"""
from fractions import Fraction as Q
import json


def check_credits(kernel_path,credit_path,j,coefficients):
    kernel=json.loads(kernel_path.read_text())
    records=json.loads(credit_path.read_text())[str(j)] if credit_path else []
    def integer(value):
        assert type(value) in (str,int)
        q=Q(value)
        assert q.denominator==1
        return q.numerator
    denominator=integer(kernel['common_denominator'])
    assert denominator>0
    credits={}
    for record in records:
        d=record['degree']
        assert type(d) is int and 1<=d<len(coefficients) and d not in credits
        block=kernel['blocks'][d]
        m=[[Q(0) for _ in range(49)] for _ in range(49)]
        for ray in block['rays']:
            weight=integer(ray['weight'])
            assert weight>0
            vector=list(map(integer,ray['vector']))
            assert len(vector)==49
            for i in range(49):
                if not vector[i]:continue
                for k in range(49):
                    if vector[k]:m[i][k]+=Q(weight,denominator)*vector[i]*vector[k]
        w=list(map(Q,record['solution']))
        assert len(w)==49
        v=[Q(int(i<7*j)) for i in range(49)]
        assert all(sum(m[i][k]*w[k] for k in range(49))==v[i] for i in range(49))
        norm=sum(x*y for x,y in zip(v,w))
        beta=Q(record['beta'])
        assert norm>0 and beta>0 and beta*norm==1
        credits[d]=beta
    assert all(c+credits.get(d,Q(0))>=0 for d,c in enumerate(coefficients) if d>=1)
    return dict(credits=credits,range_equations_checked=49*len(credits))
