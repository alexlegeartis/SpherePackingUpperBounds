#!/usr/bin/env python3
"""Replay retained Sturm proofs with Python integers only (no FLINT/SymPy).

Reconstruction from the frozen inputs is in verify.py. This independent replay
checks the certificate hashes, rational input polynomial, derivative, every
signed remainder identity, endpoint signs, zero root counts and margin claim.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import hashlib
import json
import math
import sys
from fractions import Fraction
from functools import reduce
from pathlib import Path
from time import monotonic
sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent

def trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p

def primitive(p):
    p = trim(list(p))
    if not p:
        return []
    gcd = reduce(math.gcd, map(abs,p))
    assert gcd > 0
    return [a//gcd for a in p]

def rational_primitive(p):
    den = math.lcm(*(a.denominator for a in p))
    return primitive([a.numerator*(den//a.denominator) for a in p])

def positive_pseudorem(a,b):
    # Each elimination uses |lc(b)|, never a negative multiplier. Positive
    # content removal also preserves the sign of the true rational remainder.
    r = list(a)
    while r and len(r) >= len(b):
        lc = r[-1]
        offset = len(r)-len(b)
        sign_b = 1 if b[-1] > 0 else -1
        r = [abs(b[-1])*v for v in r]
        for j,v in enumerate(b):
            r[offset+j] -= sign_b*lc*v
        r = primitive(r)
    return r

def sign_at(p,x):
    a,b=x.numerator,x.denominator
    value=p[-1]
    denominator_power=b
    for coeff in reversed(p[:-1]):
        value = value*a + coeff*denominator_power
        denominator_power *= b
    return (value>0)-(value<0)

def variations(signs):
    signs = [s for s in signs if s]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def main():
    start=monotonic()
    result=json.loads((HERE/'result.json').read_text())
    for name,digest in result['input_sha256'].items():
        raw=(HERE.parent/'density_kissing16_reweight_frozen_20261111'/name).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==digest
    raw=(HERE/'pair_coefficients.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==result['pair_coefficients_sha256']
    pair=list(map(Fraction,json.loads(raw)['coefficients']))
    left,right=map(Fraction,result['interval_exact'])
    rows=[]
    for key,shift in [('pair_sturm',Fraction(0)),('pair_with_margin_sturm',Fraction(9,10000000))]:
        proof=result[key]
        raw=(HERE/proof['chain_file']).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==proof['chain_sha256']
        chain=[list(map(int,q)) for q in json.loads(raw)]
        p=pair.copy(); p[0]+=shift
        assert chain[0]==rational_primitive(p)
        assert chain[1]==primitive([j*a for j,a in enumerate(chain[0])][1:])
        for j in range(2,len(chain)):
            rem=positive_pseudorem(chain[j-2],chain[j-1])
            assert chain[j]==primitive([-a for a in rem]), j
        assert len(chain[-1])==1 and chain[-1][0]!=0
        assert all(len(a)>len(b) for a,b in zip(chain,chain[1:]))
        signs_left=[sign_at(p,left) for p in chain]
        signs_right=[sign_at(p,right) for p in chain]
        assert signs_left==proof['left_signs'] and signs_right==proof['right_signs']
        assert signs_left[0]==signs_right[0]==-1
        va,vb=variations(signs_left),variations(signs_right)
        assert va==vb==17
        assert [len(p)-1 for p in chain]==proof['chain_degrees']
        rows.append({'certificate':key,'left_variations':va,'right_variations':vb,'roots':va-vb,'all_remainder_identities_verified':True})
    # Verify reported decimal enclosure against the exact objective.
    objective=Fraction(result['objective_exact'])
    enc=result['objective_enclosure']
    assert Fraction(enc['lower'])<=objective<=Fraction(enc['upper'])
    assert Fraction(result['strict_uniform_pair_upper_bound_exact'])==-Fraction(9,10000000)
    out={'status':'PASS','arithmetic':'Python standard-library integers and Fraction; no polynomial libraries','certificates':rows,'runtime_seconds':monotonic()-start}
    (HERE/'replay_result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':
    main()
