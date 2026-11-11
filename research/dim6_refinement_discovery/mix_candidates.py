"""Form a rational convex combination; requires fresh full verification."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('left', type=Path)
ap.add_argument('right', type=Path)
ap.add_argument('weight', type=F)
ap.add_argument('output', type=Path)
args = ap.parse_args()
a = json.loads(args.left.read_text())[0]
b = json.loads(args.right.read_text())[0]
t = args.weight
assert 0 <= t <= 1
for key in ['dimension', 'radial_degree', 'radial_nodes', 'angular_degree',
            'exclusion_squared', 'truncation_squared', 'count_squared_radii',
            'count_bounds']:
    assert a[key] == b[key], key

def mix(x, y):
    if isinstance(x, list):
        assert len(x) == len(y)
        return [mix(xx, yy) for xx, yy in zip(x, y)]
    return str((1 - t) * F(str(x)) + t * F(str(y)))

for key in ['coefficients', 'matrices', 'count_weights']:
    a[key] = mix(a[key], b[key])
a['certified'] = False
a['warning'] = 'Rational convex combination; complete fresh verification required.'
a['mixture'] = {'left': args.left.name, 'right': args.right.name,
                'right_weight': str(t)}
for key in ['objective', 'density_candidate', 'max_constraint_violation',
            'status', 'history', 'solver_iterations', 'solver_seconds']:
    a.pop(key, None)
a['separation_cuts'] = []
args.output.write_text(json.dumps([a], indent=2) + '\n')
