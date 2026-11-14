"""Screen every distinct one-row Gram pattern in a finite P48n vector pool.

Normal-row enumeration is finite and complete for the specified seed(s)/pool.
The finite dual-cluster search is heuristic; no maximum-clique claim is made.
Only exact, primitive, strictly improving witnesses are written as candidates.
"""
from pathlib import Path
from fractions import Fraction
import argparse, json, math, re, time
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from fpylll import Enumeration, EnumerationError, GSO, IntegerMatrix

P = Path(__file__).resolve().parent
ROOT = P.parents[2]
ap = argparse.ArgumentParser()
ap.add_argument('--seed-limit', type=int, default=1)
ap.add_argument('--pool', default='P48n_extended_normal_pool.npz')
ap.add_argument('--tag', default='antipode_one_row_43')
args = ap.parse_args()
t0 = time.time()
j = np.load(P / args.pool)
G = j['G'].astype(np.int64)
V = j['V'].astype(np.int16)
norms = j['norms'].astype(np.int16)
VG = V.astype(np.int64) @ G
assert np.abs(VG).max() < 32768
VG = VG.astype(np.int16)
target2 = Fraction(137**43, 2**223 * 3**42)
effective_det = float(Fraction(3, 2)**43 / target2)
EXT = ROOT / 'research/density_breakthroughs_20261113/seed_deformations/external_sphere_packing/packings_out/records/dim43-delta226'

def matrix_after(text, key):
    p = text.index(key)
    left = text.index('[[', p)
    right = text.index(']]', left) + 2
    return np.array([[int(z) for z in re.findall(r'-?\d+', row)]
                     for row in text[left:right].splitlines()], dtype=np.int64)

def evaluate(K):
    det = int(round(np.linalg.det(K)))
    if det <= 0:
        return 0.0, None
    J = np.rint(np.linalg.inv(K) * det).astype(np.int64)
    if not np.array_equal(J @ K, det * np.eye(5, dtype=np.int64)):
        return 0.0, None
    baseline = math.sqrt(effective_det / det)
    M = GSO.Mat(IntegerMatrix.from_matrix(J.tolist()), gram=True, update=True)
    try:
        shell = Enumeration(M, nr_solutions=128).enumerate(
            0, 5, min(float(np.diag(J).min()) * 2, det * .9), 0)
    except EnumerationError:
        return baseline, (det, J, np.zeros((1, 5), dtype=np.int64), 0)
    W = np.array([np.rint(row[1]).astype(np.int64) for row in shell])
    W = np.unique(np.vstack([np.zeros((1, 5), dtype=np.int64), W, -W]), axis=0)
    zero = int(np.flatnonzero(np.all(W == 0, axis=1))[0])
    A = W @ J @ W.T
    diagonal = np.diag(A)
    D = diagonal[:, None] + diagonal[None, :] - 2*A
    best = baseline
    bestW = W[[zero]]
    bestdiam = 0
    for seed in np.argsort(diagonal)[:min(36, len(W))]:
        chosen = [zero] if seed == zero else [zero, int(seed)]
        used = np.zeros(len(W), dtype=bool)
        used[chosen] = True
        cost = D[chosen].max(axis=0)
        diameter = int(cost[chosen].max())
        for size in range(len(chosen), min(18, len(W)) + 1):
            if size > len(chosen):
                options = cost.copy()
                options[used] = np.iinfo(np.int64).max
                q = int(np.argmin(options))
                chosen.append(q)
                used[q] = True
                diameter = max(diameter, int(cost[q]))
                cost = np.maximum(cost, D[q])
            beta = diameter/det
            if beta >= 6:
                break
            score = size * (1-beta/6)**21.5 * baseline
            if score > best + 1e-10:
                best, bestW, bestdiam = score, W[chosen].copy(), diameter
    return best, (det, J, bestW, bestdiam)

def check_and_save(B, K, data, tag):
    det, J, W, diameter = data
    assert int(sp.Matrix(K).det()) == det
    snf = smith_normal_form(sp.Matrix(B), domain=sp.ZZ)
    index = math.prod(abs(int(snf[i, i])) for i in range(5))
    if index != 1:
        return None
    beta = Fraction(max(int((x-y) @ J @ (x-y)) for x in W for y in W), det)
    ratio2 = Fraction(len(W)**2, det) * ((6-beta)/4)**43 / target2
    if ratio2 <= 1:
        return None
    row = dict(dimension=43, normal_determinant=det,
               saturation_index=index, cluster_size=len(W),
               cluster_diameter_squared=str(beta),
               density_ratio_squared=str(ratio2),
               density_ratio=math.sqrt(float(ratio2)),
               center_density=math.sqrt(float(target2*ratio2)),
               minimum_dependency='Published minimum 6 of P48n; not re-enumerated here.')
    np.savez(P/f'{tag}.npz', G=G, B=B, K=K, dual_integer_coordinates=W, adjugate=J)
    (P/f'{tag}.json').write_text(json.dumps(row, indent=2))
    return row

files = sorted(EXT.glob('P48n-*.cert.txt'))[:args.seed_limit]
cache = {}
groups = []
candidates = []
best = 1.0
for seednumber, file in enumerate(files):
    text = file.read_text()
    assert np.array_equal(matrix_after(text, 'ambient Gram G'), G)
    B = matrix_after(text, 'coefficient rows w.r.t.')
    K0 = B @ G @ B.T
    control, _ = evaluate(K0)
    assert abs(control - 1) < 1e-8
    print('CONTROL', file.name, control, flush=True)
    for omitted in range(5):
        keep = [i for i in range(5) if i != omitted]
        C = B[keep].astype(np.int16)
        assert int(np.abs(VG).max()) * int(np.abs(C).sum(axis=1).max()) < 32768
        dots = VG @ C.T
        # Opposite normal rows produce isometric Gram matrices.
        first = np.argmax(dots != 0, axis=1)
        signs = np.sign(dots[np.arange(len(V)), first])
        signs[signs == 0] = 1
        dots *= signs[:, None]
        signatures = np.column_stack([norms, dots])
        patterns, representatives = np.unique(signatures, axis=0, return_index=True)
        K4 = C.astype(np.int64) @ G @ C.T
        inverse4 = np.linalg.inv(K4)
        residual = patterns[:, 0] - np.einsum('bi,ij,bj->b', patterns[:, 1:], inverse4, patterns[:, 1:])
        good = residual > 1e-8
        patterns, representatives, residual = patterns[good], representatives[good], residual[good]
        order = np.argsort(residual)
        evaluated = 0
        for position in order:
            pattern = patterns[position]
            representative = int(representatives[position])
            K = np.empty((5, 5), dtype=np.int64)
            K[:4, :4] = K4
            K[4, 4] = pattern[0]
            K[:4, 4] = pattern[1:]
            K[4, :4] = pattern[1:]
            key = K.tobytes()
            if key in cache:
                continue
            score, data = evaluate(K)
            cache[key] = score
            evaluated += 1
            if score > best + 1e-8 and data is not None:
                Bnew = np.vstack([C, signs[representative] * V[representative]]).astype(np.int64)
                assert np.array_equal(Bnew @ G @ Bnew.T, K)
                row = check_and_save(Bnew, K, data, f'{args.tag}_candidate')
                if row is not None:
                    candidates.append(row)
                    best = row['density_ratio']
                    print('CANDIDATE', row, flush=True)
        group = dict(seed=file.name, omitted_row=omitted, distinct_patterns=len(patterns),
                     new_gram_evaluations=evaluated, best_ratio=best,
                     seconds=time.time()-t0)
        groups.append(group)
        print('GROUP', json.dumps(group), flush=True)
summary = dict(dimension=43, reference_center_density_squared=str(target2),
               pool=args.pool, pool_size=len(V), seed_count=len(files),
               distinct_gram_evaluations=len(cache), best_ratio=best,
               best_discovery_score=max(cache.values()), new_verified_candidates=candidates,
               groups=groups, seconds=time.time()-t0,
               scope='All distinct one-row Gram patterns for these seeds in this finite pool; dual enumeration capped at 128 rays and cluster search greedy, up to 18 points.')
(P/f'{args.tag}_summary.json').write_text(json.dumps(summary, indent=2))
print('DONE', json.dumps(summary), flush=True)
