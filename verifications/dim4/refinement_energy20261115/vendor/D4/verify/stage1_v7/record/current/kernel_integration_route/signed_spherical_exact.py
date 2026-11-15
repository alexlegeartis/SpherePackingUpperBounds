"""Exact transfer and polynomial integration of proved spherical local energies.

For a proved all-cardinality spherical energy B*N-sum_ij g(t_ij)>=0,
the Euclidean local energy has origin0, h=B*chi/2, and K=-g*chi*chi.
Its pair increment is kappa*chi, since B-g(1)=kappa. The signed cyclic
term is inserted into the polynomial numerator; it is never discarded as
if it were a nonnegative correction.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import itertools
import json
import math
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scalar_from_certificate(data, gram_blocks):
    assert data.get('verification_goal') == 'local_energy'
    assert not data.get('fixed_count'), 'Fixed-count positivity cannot transfer to arbitrary neighborhoods.'
    scale = int(data['scalar_scale'])
    assert scale > 0
    g = [Q(int(a), scale) for a in data['scalar']]
    decrease = Q(0)
    if 'uniform_repair' in data:
        repair = data['uniform_repair']
        q = list(map(Q, repair['gegenbauer_coefficients']))
        epsilon, eta = Q(repair['scale']), Q(repair['pair_shift'])
        assert q and q[0] > 0 and min(q) >= 0 and epsilon > 0 and eta >= 0
        decrease = epsilon*(sum(q)/q[0]-1)+eta
        g[0] -= decrease
    feature_scale = int(data['feature_scale'])
    assert feature_scale > 0
    kappa = Q(int(gram_blocks(data)[0][0, 0]), feature_scale**2)
    assert kappa >= 0
    return g, kappa, sum(g)+kappa, decrease


def validate_registry(records, core, folder, spherical, replay=True):
    """Recompute canonical fields, bind bytes, and freshly replay local proofs."""
    assert isinstance(records, list)
    source_folder = Path(folder)/'spherical_energy_sources'
    replay_folder = Path(folder)/'spherical_energy_replay'
    if replay and records:
        replay_folder.mkdir(exist_ok=True)
    checked, reports, dependencies = [], [], []
    seen = set()
    for record in records:
        assert record['status'] == 'EXACT_SPHERICAL_LOCAL_ENERGY_RAY'
        beta = Q(record['weight'])
        R, s = Q(record['R']), Q(record['cosine_cutoff'])
        assert beta >= 0 and 1 < R < core.BREAKS[-1] and R in core.BREAKS
        assert R*R-R-1 < 0 and 1-Q(1)/(2*R*R) <= s
        name = record['source_certificate']
        assert Path(name).name == name and name not in seen
        seen.add(name)
        path = source_folder/name
        digest = sha(path)
        assert digest == record['source_sha256']
        data = json.loads(path.read_text())
        assert Q(data['s']) == s
        g, kappa, B, decrease = scalar_from_certificate(data, spherical.gram_blocks)
        assert record['scalar_basis'] == 'U_j(t)/(j+1)'
        assert list(map(Q, record['scalar_coefficients'])) == g
        assert list(map(Q, record['global_kernel_angular_coefficients'])) == [-a for a in g]
        assert Q(record['kappa_s']) == Q(record['global_pair_cost']) == kappa
        assert Q(record['B']) == B
        assert Q(record['global_origin_coupling_coefficient']) == B/2
        assert Q(record['global_origin']) == 0
        assert Q(record['scalar_constant_decrease']) == decrease
        assert B-sum(g) == kappa
        assert len(g)-1 <= core.gram()[2], 'Homogenization degree must include the signed polynomial.'
        report = None
        if replay:
            target = replay_folder/name
            target.write_bytes(path.read_bytes())
            report = spherical.verify(target)
            assert report['status'] == 'EXACT_SPHERICAL_LOCAL_ENERGY_VERIFIED'
            assert report['certificate_sha256'] == digest
            assert Q(report['B']) == B and Q(report['cosine_cutoff']) == s
        flags = [int(core.BREAKS[j+1] <= R) for j in range(core.NS)]
        assert all(core.BREAKS[j+1] <= R or core.BREAKS[j] >= R for j in range(core.NS))
        checked.append(dict(weight=beta, R=R, s=s, g=g, B=B, kappa=kappa,
                            flags=flags, source_sha256=digest))
        reports.append(dict(source_certificate=name, source_sha256=digest,
                            weight=str(beta), global_origin='0', pair_cost=str(beta*kappa),
                            source_verification=report,
                            status='FRESH_EXACT_LOCAL_ENERGY_REPLAY_PASSED' if replay else
                                   'SCHEMA_AND_TRANSFER_IDENTITIES_ONLY_NOT_A_PROOF'))
        dependencies.append((path, digest))
    return checked, reports, dependencies


def power_coefficients(gegenbauer):
    """Normalized S3 polynomial via its closed integer coefficient formula."""
    result = [Q(0)]*len(gegenbauer)
    for ell, coefficient in enumerate(gegenbauer):
        for j in range(ell//2+1):
            k = ell-2*j
            result[k] += coefficient*Q((-1)**j*math.comb(ell-j, j)*2**k, ell+1)
    return result


def direct_g(g, t):
    """Independent three-term recurrence, including the diagonal t=1."""
    p0, result = Q(1), g[0]
    if len(g) == 1:
        return result
    p1 = t
    result += g[1]*p1
    for ell in range(1, len(g)-1):
        p2 = ((2*ell+2)*t*p1-ell*p0)/(ell+2)
        result += g[ell+1]*p2
        p0, p1 = p1, p2
    return result


def direct_cyclic(kind, point, regions, rays, cutoff):
    a, b, c = map(Q, point)
    result = Q(0)
    if kind == 'full':
        sides = [a, b, c]
        pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
        for ray in rays:
            for i, j, k in pairs:
                if ray['flags'][regions[i]] and ray['flags'][regions[j]]:
                    u, v, w = sides[i], sides[j], sides[k]
                    result -= ray['weight']*direct_g(ray['g'], (u*u+v*v-w*w)/(2*u*v))
    else:
        assert kind == 'outside'
        t = -1+c*((a+b)**2-cutoff**2)/(2*a*b)
        for ray in rays:
            if ray['flags'][regions[0]] and ray['flags'][regions[1]]:
                result -= ray['weight']*direct_g(ray['g'], t)
    return result


def sparse_numerator(kind, regions, rays, L, cutoff):
    """Return exactly factor*Gamma_spherical in the global branch coordinates.

    Full factor=(2abc)^L. Exterior factor=(2ab)^L and the coordinate z in
    [0,1] gives t=-1+z*((a+b)^2-cutoff^2)/(2ab). The coefficients are
    rational, and the radial flags are constant on each closed branch.
    Both one-sided branch values are retained at a support boundary.
    """
    out = {}

    def add(index, value):
        out[index] = out.get(index, Q(0))+value

    for ray in rays:
        assert ray['weight'] >= 0 and len(ray['g'])-1 <= L
        p = power_coefficients(ray['g'])
        if kind == 'full':
            for ai, aj, ak in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
                if not (ray['flags'][regions[ai]] and ray['flags'][regions[aj]]):
                    continue
                for k, pk in enumerate(p):
                    if pk == 0:
                        continue
                    pre = -ray['weight']*pk*2**(L-k)
                    for u in range(k+1):
                        for v in range(k-u+1):
                            w = k-u-v
                            index = [0, 0, 0]
                            index[ai], index[aj], index[ak] = L-k+2*u, L-k+2*v, L+2*w
                            add(tuple(index), pre*math.comb(k, u)*math.comb(k-u, v)*(-1)**w)
        else:
            assert kind == 'outside'
            if not (ray['flags'][regions[0]] and ray['flags'][regions[1]]):
                continue
            # W=-2ab+z*(a^2+2ab+b^2-cutoff^2), so t=W/(2ab).
            terms = [((1, 1, 0), Q(-2)), ((2, 0, 1), Q(1)),
                     ((1, 1, 1), Q(2)), ((0, 2, 1), Q(1)),
                     ((0, 0, 1), -cutoff**2)]
            power = {(0, 0, 0): Q(1)}
            for k, pk in enumerate(p):
                if pk:
                    pre = -ray['weight']*pk*2**(L-k)
                    for (u, v, w), coefficient in power.items():
                        add((L-k+u, L-k+v, w), pre*coefficient)
                if k+1 < len(p):
                    following = {}
                    for old, coefficient in power.items():
                        for shift, multiplier in terms:
                            index = tuple(a+b for a, b in zip(old, shift))
                            following[index] = following.get(index, Q(0))+coefficient*multiplier
                    power = {index: value for index, value in following.items() if value}
    return {index: value for index, value in out.items() if value}


def audit_branch(core, kind, regions, rays, sparse=None):
    L = core.gram()[2]
    sparse = sparse_numerator(kind, regions, rays, L, core.BREAKS[-1]) if sparse is None else sparse
    domain = [tuple(core.BREAKS[j:j+2]) for j in regions]
    if kind == 'outside':
        domain.append((Q(0), Q(1)))
    points = [(Q(1, 3), Q(2, 5), Q(3, 7)),
              (Q(0), Q(1), Q(1, 2)), (Q(1), Q(0), Q(1)),
              (Q(0), Q(0), Q(0)), (Q(1), Q(1), Q(1))]
    for unit in points:
        point = [a+(b-a)*u for (a, b), u in zip(domain, unit)]
        power_value = sum((value*math.prod(point[j]**index[j] for j in range(3))
                           for index, value in sparse.items()), Q(0))
        factor = (2*point[0]*point[1]*(point[2] if kind == 'full' else 1))**L
        expected = factor*direct_cyclic(kind, point, regions, rays, core.BREAKS[-1])
        assert power_value == expected, ('Signed polynomial identity', kind, regions, unit)
    return dict(status='EXACT_SIGNED_SPHERICAL_BRANCH_IDENTITIES_VERIFIED',
                kind=kind, regions=list(regions), checks=len(points), sparse_terms=len(sparse))


def build_signed(core, kind, regions, eta, rays):
    C, denominator, domain = core.build(kind, regions, eta=eta)
    L = core.gram()[2]
    sparse = sparse_numerator(kind, regions, rays, L, core.BREAKS[-1])
    report = audit_branch(core, kind, regions, rays, sparse)
    dd = math.lcm(denominator, *(value.denominator for value in sparse.values()))
    answer = C.copy()*(dd//denominator)
    for index, value in sparse.items():
        assert all(0 <= k < n for k, n in zip(index, answer.shape))
        answer[index] += int(value*dd)
    return answer, dd, domain, report
