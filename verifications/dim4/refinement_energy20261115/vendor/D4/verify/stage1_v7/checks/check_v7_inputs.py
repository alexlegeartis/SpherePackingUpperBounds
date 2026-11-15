#!/usr/bin/env python3
"""Check the exact V7 data and replay the archived exact scalar Sturm code.

Uses the Python standard library; no archived success flag is a scalar sign
test. The positivity of every Gram weight gives a constructive PSD proof.
This checker does not verify global pair or cyclic inequalities.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import importlib.util
import json
import math
import sys

sys.dont_write_bytecode = True
CANDIDATE = Path('current/candidates/round9_six_graphs_degree18')
SCALAR_VERIFIER_SHA256 = '069596bffc7444e6f739cfadf9c70f4ddc4a650f0a2632622e904e0dde2d61e7'


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(record):
    record = record.resolve()
    source = record / CANDIDATE
    kernel = json.loads((source/'kernel.json').read_text())
    origin = json.loads((source/'origin.json').read_text())
    scalar = json.loads((source/'scalar.json').read_text())
    graph = json.loads((source/'graph_energies.json').read_text())
    require(kernel['dimension'] == scalar['dimension'] == 4, 'Dimension mismatch')
    require(kernel['local_degree'] == 6 and kernel['angular_degree'] == 18, 'Degrees differ')
    edges = list(map(Q, kernel['breaks']))
    bases = [list(map(Q, row)) for row in kernel['radial_basis_intervals']]
    require(len(edges) == 8 and edges[0] == 1 and edges[-1] == Q(7,4), 'Radial support differs')
    require(len(bases) == 7 and len(kernel['blocks']) == 19, 'Incomplete basis or Gram data')
    for i, (left, right) in enumerate(bases):
        require(left <= edges[i] < edges[i+1] <= right, 'Invalid radial mother interval')
    den = int(kernel['common_denominator'])
    require(den > 0, 'Nonpositive Gram denominator')
    counts = []
    for ell, block in enumerate(kernel['blocks']):
        n = 50 if ell == 0 else 49
        require(block['dimension'] == n, 'Wrong Gram block dimension')
        for ray in block['rays']:
            require(type(ray['weight']) is int and ray['weight'] > 0, 'Invalid Gram ray weight')
            require(len(ray['vector']) == n and all(type(v) is int for v in ray['vector']), 'Invalid Gram ray vector')
        counts.append(len(block['rays']))
    require(sum(counts) == 628, 'Unexpected ray count')
    v, u = (list(map(Q, origin[name])) for name in ('v', 'u'))
    require(len(v) == len(u) == len(origin['pair_origin_v']) == 49, 'Missing origin coefficients')
    require(len(origin['functional_origin_v']) == 15, 'Missing physical functional coefficients')
    # Compute M*u directly from its constructive sum of outer products.
    product = [Q(0) for _ in range(49)]
    for ray in kernel['blocks'][0]['rays']:
        q = ray['vector'][1:]
        factor = Q(ray['weight'], den) * sum((a*b for a,b in zip(q,u)), Q(0))
        for i, a in enumerate(q):
            product[i] += a*factor
    require(product == v, 'Mu=v fails')
    require(sum((a*b for a,b in zip(v,u)), Q(0)) == Q(origin['kappa']) >= 0, 'Completion constant fails')
    f, g = (list(map(Q, scalar[key])) for key in ('F','G'))
    require(len(f) == len(g) == 50 and Q(scalar['alpha']) == 2, 'Scalar degree or scale differs')
    require(Q(scalar['cutoff']) == Q(7,4), 'Scalar cutoff differs')
    transformed = [sum((f[k]*(-1)**j*math.factorial(k)*math.comb(k+1,k-j)/math.factorial(j)
                        for k in range(j,len(f))), Q(0)) for j in range(len(f))]
    require(transformed == g and g[0] == 1, 'Fourier polynomial identity fails')
    # Bind the reused exact scalar algorithm to the preserved source inventory.
    scalar_code = source/'verify_scalar.py'
    # Stage 1 export: the expected digest is the entry for this file in the
    # original record/MANIFEST.json, pinned here because that 2,957-file
    # manifest is not exported.
    require(sha(scalar_code) == SCALAR_VERIFIER_SHA256, 'Scalar verifier changed')
    spec = importlib.util.spec_from_file_location('v7_scalar_exact_replay', scalar_code)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    signs = {}
    for label, poly, left, wanted in [('spatial_tail', f, Q(49,8), -1), ('Fourier_half_line', g, Q(0), 1)]:
        sequence = module.sturm(module.integers(poly))
        a = [module.sign_at(p,left) for p in sequence]
        b = [module.sign_at(p,None) for p in sequence]
        va, vb = module.variations(a), module.variations(b)
        require(va == vb and a[0] == b[0] == wanted, 'Scalar half-line sign fails')
        signs[label] = {'left': str(left), 'left_variations': va, 'infinity_variations': vb,
                        'root_count': va-vb, 'sign': wanted,
                        'signs_at_left': a, 'signs_at_infinity': b}
    require(graph['schema'] == 'INDEPENDENT_EXTERIOR_GRAPH_ENERGIES_V7', 'Not V7')
    require(graph['base_input_sha256'] == {n: sha(source/n) for n in ('kernel.json','origin.json','scalar.json')}, 'Graph/base binding fails')
    require(not (source/'gaussian_tail_repair.json').exists(), 'Unexpected Gaussian repair')
    local_fields = ['pointed_caps','conditional_radial_rays','pointed_tents',
                    'conditional_count_polynomials','independent_antipodal_allocations',
                    'appended_spherical_local_energy_rays','appended_broad19_nonorigin_psd_rays']
    return {'status':'V7_EXACT_DATA_GRAM_ORIGIN_AND_SCALAR_CHECKS_PASSED',
            'scope':'Constructive Gram PSD, original origin equations, and complete scalar half-line signs; not a global cyclic replay.',
            'gram_ray_count':sum(counts), 'gram_rays_by_angular_degree':counts,
            'local_energy_record_counts':{key:len(origin[key]) for key in local_fields},
            'radial_breaks':list(map(str,edges)), 'radial_basis_intervals':kernel['radial_basis_intervals'],
            'physical_pair_origin_entries':49, 'functional_origin_entries':15,
            'original_origin_range_equations':49, 'scalar_coefficients_in_each_polynomial':50,
            'scalar_sturm_code_sha256':sha(scalar_code), 'scalar_half_line_proofs':signs,
            'gaussian_tail_repair_used':False,
            'input_sha256':{n:sha(source/n) for n in ('kernel.json','origin.json','scalar.json','graph_energies.json')}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output:
        args.output = args.output.resolve()
        if args.output.exists() or args.output.is_relative_to(args.record.resolve()):
            parser.error('--output must be a new file outside the preserved record')
    data = json.dumps(check(args.record), indent=2)+'\n'
    if args.output:
        args.output.write_text(data)
    print(data, end='')
