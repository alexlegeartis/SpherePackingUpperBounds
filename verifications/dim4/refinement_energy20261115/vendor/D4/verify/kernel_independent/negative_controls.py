#!/usr/bin/env python3
"""Negative controls for run.py. Each must FAIL with an exact mismatch.

    python3 negative_controls.py --out DIR

Three runs on the same small case subset (the pair pieces, one full case and
one exterior case), each in its own subfolder of DIR:
  baseline        unchanged; every compared object must be equal
  implementation  kernel.py with one signed spherical weight times (1 + 1e-9)
  data            kernel.py reading a copy of data/v7 in which the weight of
                  ray 0 of angular block 3 in kernel.json is increased by 1
                  (stage 1 still reads data/v7)
A control passes when its run reports an exact polynomial mismatch in the
pair term or in a compared case while the baseline reports equality there.
DIR/CONTROLS.json records the outcome; exit 0 if both controls behave.
"""
from pathlib import Path
import argparse
import copy
import json
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]
SUBSET = 'full:0,0,0;outside:5,5;full:3,5,6'


def tensor_input_controls():
    """Regression: malformed captured coefficients must not be truncated."""
    from fractions import Fraction as Q
    import flint
    import numpy as np
    from extract_stage1 import compare_polynomial
    ctx = flint.fmpq_mpoly_ctx.get(('x0','x1','x2'), 'lex')
    one = ctx.from_dict({(0,0,0):1})
    assert compare_polynomial(np.array([[[1]]],dtype=object),1,one)['equal']
    cases = [
        ('fractional_coefficient',np.array([[[Q(3,2)]]],dtype=object),1),
        ('fractional_denominator',np.array([[[1]]],dtype=object),Q(3,2)),
        ('zero_denominator',np.array([[[1]]],dtype=object),0),
        ('wrong_tensor_rank',np.array([[1]],dtype=object),1),
    ]
    result = {}
    for name, coefficients, denominator in cases:
        try:
            compare_polynomial(coefficients,denominator,one)
        except ValueError as error:
            result[name] = dict(status='REJECTED', reason=str(error))
        else:
            result[name] = dict(status='NOT_REJECTED')
    return result


def invalid_data_controls(out):
    """Each malformed consumer input must be rejected before coercion."""
    import kernel as km
    data = {f.name: json.loads(f.read_text()) for f in (PKG/'data'/'v7').glob('*.json')}
    mutations = [
        ('fractional_degree', 'kernel.json', lambda x: x.update(local_degree='13/2')),
        ('fractional_denominator', 'kernel.json', lambda x: x.update(common_denominator='3/2')),
        ('boolean_dimension', 'kernel.json', lambda x: x.update(dimension=True)),
        ('fractional_ray', 'kernel.json', lambda x: x['blocks'][1]['rays'][0]['vector'].__setitem__(0,'1/2')),
        ('short_origin_row', 'origin.json', lambda x: x['v'].pop()),
        ('fractional_cap_count', 'origin.json', lambda x: x['pointed_caps'][0].update(M='3/2')),
        ('negative_cap_weight', 'origin.json', lambda x: x['pointed_caps'][0].update(weight='-1')),
        ('duplicate_count_radius', 'origin.json', lambda x: x['nested_count_quadratic']['radii'].__setitem__(1,x['nested_count_quadratic']['radii'][0])),
        ('fractional_gauge_count', 'origin.json', lambda x: x['full_pair_count_gauge'].update(count_bound='49/2')),
        ('float_scalar', 'scalar.json', lambda x: x.update(alpha=2.0)),
    ]
    destination = out/'invalid_data'
    destination.mkdir()
    results = {}
    for name, filename, mutate in mutations:
        changed = copy.deepcopy(data)
        mutate(changed[filename])
        for f, value in changed.items():
            (destination/f).write_text(json.dumps(value))
        try:
            km.Kernel(destination)
        except (AssertionError, ValueError, TypeError) as error:
            results[name] = dict(status='REJECTED', reason=type(error).__name__+': '+str(error))
        else:
            results[name] = dict(status='NOT_REJECTED')
    return results


def run(out, extra):
    argv = [sys.executable, str(HERE/'run.py'), '--out', str(out), '--only', SUBSET]+extra
    code = subprocess.run(argv, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT).returncode
    result = json.loads((out/'RESULT.json').read_text())
    compared = [c for c in result.get('cases', []) if c['status'] != 'MISSING']
    return dict(exit=code, status=result['status'],
                pair=(result.get('pair') or {}).get('status'),
                cases={'%s %s' % (c['kind'], c['regions']): c['status'] for c in compared},
                first_mismatch=next((c['failure'] for c in compared if c['status'] != 'PASS'), None))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {}
    report['invalid_inputs'] = invalid_data_controls(out)
    report['invalid_tensor_inputs'] = tensor_input_controls()
    report['baseline'] = run(out/'baseline', [])
    report['implementation'] = run(out/'implementation', ['--perturb', 'implementation'])
    data = out/'data_v7_perturbed'
    shutil.copytree(PKG/'data'/'v7', data)
    kernel = json.loads((data/'kernel.json').read_text())
    kernel['blocks'][3]['rays'][0]['weight'] += 1
    (data/'kernel.json').write_text(json.dumps(kernel))
    report['data'] = run(out/'data', ['--kernel-data', str(data)])
    base = report['baseline']
    base_ok = base['pair'] == 'PASS' and len(base['cases']) == 3 and all(v == 'PASS' for v in base['cases'].values())
    controls = {}
    for name in ['implementation', 'data']:
        r = report[name]
        caught = r['status'] == 'FAIL' and (r['pair'] == 'FAIL' or any(v == 'FAIL' for v in r['cases'].values()))
        controls[name] = 'FAILS_AS_REQUIRED' if caught else 'NOT_CAUGHT'
    report['baseline_subset_equal'] = base_ok
    report['controls'] = controls
    malformed_ok = all(v['status'] == 'REJECTED'
                       for group in ('invalid_inputs','invalid_tensor_inputs')
                       for v in report[group].values())
    report['status'] = 'PASS' if base_ok and malformed_ok and all(v == 'FAILS_AS_REQUIRED' for v in controls.values()) else 'FAIL'
    (out/'CONTROLS.json').write_text(json.dumps(report, indent=1)+'\n')
    print(json.dumps(report, indent=1)[:3000])
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
