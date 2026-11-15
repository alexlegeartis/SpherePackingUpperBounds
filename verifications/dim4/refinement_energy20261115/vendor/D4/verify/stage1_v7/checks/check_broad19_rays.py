#!/usr/bin/env python3
"""Check that the 5 broad19 nonorigin PSD rays of origin.json are Gram rays of kernel.json.

origin.json lists five rank-one components `appended_broad19_nonorigin_psd_rays`
in a 19-coordinate basis: two degree-six Bernstein blocks on the mother
intervals [1,3/2] and [3/2,7/4], then five count indicators at the radii
R_1..R_5. No verifier reads these records; their kernel contribution must
already be inside kernel.json. This script lifts each 19-vector to the
49-coordinate piecewise basis of kernel.json (piece i uses the block of its
mother interval plus every indicator whose radius is at least the piece's
right end) and finds a ray in the same angular block whose integer vector is
a positive rational multiple of the lift. It also checks origin_value = 0.
Standard library only; exact arithmetic.

Usage: check_broad19_rays.py RECORD --output FILE
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json

CANDIDATE = 'current/candidates/round9_six_graphs_degree18'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def lift(coefficients, breaks, radii):
    require(len(coefficients) == 19 and len(radii) == 5, 'broad19 vector must have 19 entries')
    out = []
    for left, right in zip(breaks, breaks[1:]):
        block = coefficients[0:7] if right <= Q(3, 2) else coefficients[7:14]
        offset = sum((c for c, r in zip(coefficients[14:], radii) if right <= r), Q(0))
        out.extend(x+offset for x in block)
    return out


def check(record):
    folder = record/CANDIDATE
    origin = json.loads((folder/'origin.json').read_text())
    kernel = json.loads((folder/'kernel.json').read_text())
    breaks = [Q(x) for x in kernel['breaks']]
    require(len(breaks) == 8 and breaks[6] == Q(3, 2) and breaks[7] == Q(7, 4), 'unexpected breaks')
    require(kernel['local_degree'] == 6, 'unexpected radial degree')
    rows = []
    records = origin['appended_broad19_nonorigin_psd_rays']
    require(len(records) == 5, 'expected 5 broad19 rays')
    for record_ in records:
        require(record_['kind'] == 'fixed_broad19_rank_one_psd', 'unexpected kind')
        require(Q(record_['origin_value']) == 0, 'nonzero origin value')
        require([Q(x) for x in record_['indicator_radii']] == breaks[1:6], 'indicator radii differ from breaks')
        require([list(map(Q, b)) for b in record_['broad_intervals']] == [[Q(1), Q(3, 2)], [Q(3, 2), Q(7, 4)]],
                'unexpected broad intervals')
        degree = record_['angular_degree']
        require(type(degree) is int and 1 <= degree <= 18, 'angular degree out of range')
        lifted = lift([Q(x) for x in record_['coefficients']], breaks, [Q(x) for x in record_['indicator_radii']])
        require(any(lifted), 'zero ray')
        pivot = next(i for i, x in enumerate(lifted) if x != 0)
        block = kernel['blocks'][degree]
        found = None
        for index, ray in enumerate(block['rays']):
            vector = ray['vector']
            require(len(vector) == len(lifted), 'dimension mismatch')
            if vector[pivot] == 0:
                continue
            scale = Q(vector[pivot])/lifted[pivot]
            if scale > 0 and all(Q(a) == scale*b for a, b in zip(vector, lifted)):
                require(type(ray['weight']) is int and ray['weight'] > 0, 'nonpositive Gram weight')
                found = dict(ray_index=index, vector_over_lift=str(scale), gram_weight_positive=True)
                break
        require(found is not None, f'broad19 ray of angular degree {degree} is not a kernel ray')
        rows.append(dict(angular_degree=degree, **found))
    return {'status': 'BROAD19_RAYS_ARE_KERNEL_GRAM_RAYS',
            'scope': 'Each broad19 nonorigin PSD component is a positive multiple of one integer Gram ray in kernel.json; its PSD energy is part of the 628-ray Gram kernel.',
            'rays': rows}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = check(args.record.resolve())
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
