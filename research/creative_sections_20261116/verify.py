"""Exact 57D lamination, 78D/79D antipodes, and a rejected 56D code.

Reuses the source-bound MW56/MW80 reconstructions and the completed MW56
coset certificate. --export emits full period bases and rational translates.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse, json, sys
import numpy as np
from flint import fmpz_mat, fmpq_mat, nmod_mat

if not __debug__:
    raise RuntimeError('Run without -O.')
P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent/'cm58_20261116'))
import verify as cm
import verify_frontier as frontier
PIN = 'f0217194ec528125f6b115b254584b0439c68a5a'


def decsqrt(x):
    with localcontext() as ctx:
        ctx.prec = 70
        return str((Decimal(x.numerator)/x.denominator).sqrt())


def gain(x):
    with localcontext() as ctx:
        ctx.prec = 70
        g = 100*(Decimal(decsqrt(x))-1)
        low = F(int(g*10**6), 10**6)
    high = low+F(1,10**6)
    assert (1+low/100)**2 <= x < (1+high/100)**2
    return {'percent':str(g), 'interval':[str(low),str(high)],
            'display':f'{g:.6f}%'}


def audit(export=None):
    if export:
        export.mkdir(parents=True, exist_ok=False)
    checkpoint = json.loads((P/'checkpoint.json').read_text())
    assert checkpoint['commit'] == PIN
    old = {a['dimension']:a for a in checkpoint['rows']}
    assert set(old) == {57,78,79}
    result58, full58 = cm.audit()
    assert result58['enumeration']['complete_full_runs'] == 2
    # This is the leading principal 57D lattice of the certified 58D Gram.
    g57 = fmpq_mat([r[:57] for r in full58.tolist()[:57]])
    assert F(str(g57.det())) == F(7,2)*29**8
    d57 = F(5,2)**57/F(str(g57.det()))
    assert F(13,2)+F(7,2) == 10 and 4*F(7,2) >= 10
    ext57 = F(5,2)**56/F(4*29**8)
    assert F(old[57]['comparison']['density_squared']) == ext57
    rows = [{'dimension':57, 'name':'Explicit MW56 deep half-coset lamination',
             'kind':'new_explicit_lattice', 'density_squared':str(d57),
             'center_density':decsqrt(d57), 'minimum_squared':10,
             'gram_determinant':str(g57.det()),
             'comparison_name':'Shioda MW56 x minimum-matched A1 (derived)',
             'comparison_density_squared':str(ext57),
             'external_gain':gain(d57/ext57),
             'incremental_gain':gain(d57/F(old[57]['best_checked']['density_squared'])),
             'half_coset_squared_depth_at_least':'13/2',
             'vertical_squared_height':'7/2'}]
    assert g57[0,0] == 10 or any(g57[i,i] == 10 for i in range(56))
    grams = {'gram57.json':[[str(x) for x in r] for r in g57.tolist()]}
    G, _ = frontier.parent80()
    assert G.det() == 41**12
    inv = G.inv()
    source80 = json.loads((P.parent/'cross_dimension_lower_20261115/l80_results.json').read_text())
    bare = {a['dimension']:F(a['bare_section_ceiling_center_density_squared_exact'])
            for a in source80['rows']}
    for a in json.loads((P/'sections80.json').read_text())['rows']:
        n = a['dimension']; assert n in (78,79)
        k = 80-n; B = fmpz_mat(a['normal_labels']); S = fmpz_mat(a['cluster_labels'])
        assert B.nrows() == k and B.ncols() == 80 and S.ncols() == k
        assert len({tuple(r) for r in S.tolist()}) == S.nrows() == k+1
        A, T = B.transpose().hnf(transform=True)
        assert abs(T.det()) == 1
        assert A == T*B.transpose()
        assert A.tolist()[:k] == [[int(i==j) for j in range(k)] for i in range(k)]
        lifts, N = fmpz_mat(T.tolist()[:k]), fmpz_mat(T.tolist()[k:])
        assert N*B.transpose() == fmpz_mat(n,k)
        assert lifts*B.transpose() == fmpz_mat(np.eye(k,dtype=int).tolist())
        K = fmpq_mat(B)*inv*fmpq_mat(B.transpose())
        Q = N*G*N.transpose()
        assert Q.det() == G.det()*K.det()
        shifts = fmpq_mat(S*lifts)-fmpq_mat(S)*K.inv()*fmpq_mat(B)*inv
        assert shifts*fmpq_mat(B.transpose()) == fmpq_mat(S.nrows(),k)
        coords = shifts*fmpq_mat(G*N.transpose())*Q.inv()
        assert coords*fmpq_mat(N) == shifts
        beta = F(0)
        for i in range(S.nrows()):
            for j in range(i):
                v = fmpq_mat([[int(S[i,t]-S[j,t]) for t in range(k)]])
                beta = max(beta,F(str((v*K.inv()*v.transpose())[0,0])))
        assert 0 < beta < 14
        witness = fmpz_mat([a['attainment_witness']])
        i,j = a['attainment_label_pair']
        labels = fmpz_mat([[int(S[i,t]-S[j,t]) for t in range(k)]])
        assert witness*B.transpose() == labels
        assert (witness*G*witness.transpose())[0,0] == 14
        assert F(str((fmpq_mat(labels)*K.inv()*fmpq_mat(labels.transpose()))[0,0])) == beta
        d2 = S.nrows()**2*((14-beta)/4)**n/F(int(Q.det()))
        ext = {78:F(2**74,3),79:F(2**76)}[n]
        assert ext == bare[n] == F(old[n]['comparison']['density_squared'])
        assert d2 > F(old[n]['best_checked']['density_squared']) > ext
        rows.append({'dimension':n, 'name':f'Explicit MW80 antipode, {S.nrows()} translates',
                     'kind':'new_explicit_periodic_packing', 'density_squared':str(d2),
                     'center_density':decsqrt(d2), 'parent_dimension':80,
                     'parent_minimum_squared_inherited':14, 'normal_rank':k,
                     'normal_gram':[[str(x) for x in r] for r in K.tolist()],
                     'normal_determinant':str(K.det()), 'period_determinant':str(Q.det()),
                     'cosets':S.nrows(), 'cluster_squared_diameter':str(beta),
                     'minimum_squared_guarantee':str(14-beta),
                     'minimum_attained_by_exact_parent_vector':True,
                     'comparison_name':'Bachoc-Nebe L80 / minimum-vector normal (79D) or 2 A2 normal (78D), derived',
                     'comparison_density_squared':str(ext), 'external_gain':gain(d2/ext),
                     'incremental_gain':gain(d2/F(old[n]['best_checked']['density_squared'])),
                     'complete_integer_kernel_checked':True,
                     'all_projected_shifts_checked_exactly':True})
        if export:
            (export/f'dim{n}.json').write_text(json.dumps({
                'period_basis_in_parent':N.tolist(),'period_gram':Q.tolist(),
                'shifts_in_period_coordinates':[[str(x) for x in r] for r in coords.tolist()]},
                indent=2,default=int)+'\n')
    # A finite negative certificate for deleting the known 56D shell by a code.
    rejected = json.loads((P/'rejected_code56.json').read_text())
    data = np.load(P.parent/'cm58_20261116/certificate.npz',allow_pickle=False)
    t = fmpz_mat(data['T'].tolist()); C = t.inv()*fmpz_mat(data['C'].tolist())*t
    assert all(x.denominator == 1 for r in C.tolist() for x in r)
    C = np.array([[int(x) for x in r] for r in C.tolist()],dtype=np.int64)
    checks = []
    for mask in rejected['masks']:
        v = np.array([(mask>>j)&1 for j in range(56)],dtype=np.int64)
        checks.extend([v,(C@v)%2])
    checks = np.array(checks,dtype=np.int64)
    assert nmod_mat(checks.tolist(),2).rank() == 6
    g = fmpz_mat(data['G'].tolist())
    maps = np.load(P.parent/'cm58_20261116/frontier/isometries56.npz',allow_pickle=False)['maps']
    known,_ = frontier.orbit(g,data['point_coordinates'][0],maps,10)
    assert len(known) == rejected['known_rays_avoided'] == 2436
    assert np.all(np.any((known@checks.T)%2,axis=1))
    w = fmpz_mat([rejected['witness']])
    assert int((w*g*w.transpose())[0,0]) == 10
    assert np.all((np.array(rejected['witness'],dtype=np.int64)@checks.T)%2 == 0)
    assert [a['dimension'] for a in rows] == [57,78,79]
    return {'status':'verified_three_stronger_explicit_lower_bounds', 'rows':rows,
            'all_source_and_period_checks_exact':True, 'checkpoint_commit':PIN,
            'maximum_parent_or_target_dimension':80,
            'rejected_code56':{'index':64,'known_rays_deleted':2436,
                               'exact_surviving_norm10_witness_checked':True},
            'inherited_inputs':'Published MW56/MW80 minima; MW56 interval coset certificate with its stated arithmetic assumptions.'}, grams


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--write-result',action='store_true')
    ap.add_argument('--export',type=Path); args = ap.parse_args()
    result,grams = audit(args.export)
    if args.write_result:
        (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        for name,body in grams.items(): (P/name).write_text(json.dumps(body,indent=2)+'\n')
    else:
        assert json.loads((P/'results.json').read_text()) == result
        for name,body in grams.items(): assert json.loads((P/name).read_text()) == body
    for row in result['rows']:
        print(row['dimension'],row['center_density'],'prior',row['incremental_gain']['display'],
              'external',row['external_gain']['display'])
