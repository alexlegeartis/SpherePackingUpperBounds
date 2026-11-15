#!/usr/bin/env python3
"""Quick mode only: rerun the coordinator's pre-fork checks that --stage pair skips.

The V7 coordinator runs two parent checks only when it also runs cyclic
branches: the exact spherical count certificates (certify.py) and the
scalar Sturm half-line proofs (verify_scalar.py). This script loads the
same two unchanged modules from the working record and runs them on the
same inputs. It applies the same acceptance tests the coordinator applies
before it forks workers. Full mode does not use this script.

Usage: quick_parent_checks.py RECORD --output FILE
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

CANDIDATE = 'current/candidates/round9_six_graphs_degree18'
SPHERICAL = 'current/base/new_work/spherical_sdp'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def walk(value):
    if isinstance(value, dict):
        yield value
        for item in value.values():
            yield from walk(item)
    elif isinstance(value, list):
        for item in value:
            yield from walk(item)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Run without -O.')
    record = args.record.resolve()
    folder = record/CANDIDATE
    here = record/SPHERICAL
    sys.path.insert(0, str(here))
    sys.path.insert(0, str(folder))
    origin = json.loads((folder/'origin.json').read_text())

    # Every certificate the origin data cites by name, with its pinned hash.
    cited = {}
    radius_checks = []
    for item in walk(origin):
        for name_key, hash_key in [('certificate', 'certificate_sha256'),
                                   ('spherical_certificate', 'spherical_certificate_sha256')]:
            if name_key in item and hash_key in item:
                name, digest = item[name_key], item[hash_key]
                assert isinstance(name, str) and Path(name).name == name
                assert cited.setdefault(name, digest) == digest
                if name_key == 'certificate' and 'R' in item and 'M' in item:
                    radius_checks.append((name, F(item['R']), item['M']))
    assert cited, 'No spherical count certificate is cited.'

    spherical = load('spherical_nested_exact_quick', here/'certify.py')
    replay = args.output.parent/'spherical_replay_quick'
    replay.mkdir(parents=True, exist_ok=False)
    reports = {}
    for name, digest in sorted(cited.items()):
        source = here/name
        assert hashlib.sha256(source.read_bytes()).hexdigest() == digest, name
        copy = replay/name
        copy.write_bytes(source.read_bytes())
        report = spherical.verify(copy)
        assert report['certificate_sha256'] == digest, name
        reports[name] = report
    # The coordinator's tests for ordinary count dependencies.
    for name, radius, maximum in radius_checks:
        report = reports[name]
        assert 1 <= radius and radius*radius-radius-1 < 0
        assert 1-F(1)/(2*radius*radius) <= F(report['cosine_cutoff'])
        assert report['count_bound'] <= maximum

    scalar = load('verify_scalar_quick', folder/'verify_scalar.py')
    scalar_report = scalar.verify()

    result = {'status': 'QUICK_PARENT_CHECKS_PASSED',
              'scope': 'Spherical count certificates and scalar half-line signs that coordinator --stage pair skips.',
              'spherical_certificates': {name: {'sha256': cited[name], 'status': r.get('status'),
                                                'count_bound': r.get('count_bound'),
                                                'cosine_cutoff': r.get('cosine_cutoff')}
                                         for name, r in reports.items()},
              'radius_count_checks': len(radius_checks),
              'scalar': scalar_report}
    args.output.write_text(json.dumps(result, indent=2, default=str)+'\n')
    print(json.dumps({'status': result['status'], 'certificates': sorted(reports),
                      'radius_count_checks': len(radius_checks)}))


if __name__ == '__main__':
    main()
