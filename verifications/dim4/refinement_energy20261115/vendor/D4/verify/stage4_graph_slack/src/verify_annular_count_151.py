"""Fresh exact 151-point count throughout [33/28,7/4], independent of E."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
import sys
from verify_exterior_graphs import parse, harmonic_power, prove_nonnegative, max_angle


def verify(data):
    assert data['schema'] == 'S3_ANNULAR_COUNT_V1'
    radial = tuple(map(parse, data['radial_interval']))
    cutoff = parse(data['spherical_cutoff'])
    assert radial == (Q(33, 28), Q(7, 4))
    assert max_angle(radial, radial) == cutoff == Q(41, 49)
    a = list(map(parse, data['harmonic_coefficients']))
    assert a[0] > 0 and min(a[1:]) >= 0
    p = harmonic_power(a)
    sign = prove_nonnegative([-v for v in p], Q(-1), cutoff)
    cap = sum(a)/a[0]
    count = cap.numerator // cap.denominator
    assert count == data['claimed_count_bound'] == 151
    return dict(status='EXACT_ANNULAR_COUNT_151_PASS',
                global_density_bound_proved=False,
                radial_interval=data['radial_interval'],
                spherical_cutoff=str(cutoff),
                count_bound=count, f0=str(a[0]), f1=str(sum(a)),
                exact_count_upper=str(cap), whole_interval_sign_proof=sign,
                inherited_count_used=False)


if __name__ == '__main__':
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('annular_count_151_certificate.json')
    report = verify(json.loads(path.read_text()))
    report['certificate_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    report['verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(report, indent=2))
