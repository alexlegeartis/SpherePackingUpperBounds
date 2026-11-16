"""Exact source-to-packing audit for the stronger MW92 simplex antipodes.

Published parent minimum16 is inherited. No enumeration or discovery output
is used. --export emits complete integer periods and rational translates.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, importlib.util, json
import numpy as np
from flint import fmpz_mat, fmpq_mat

if not __debug__:
    raise RuntimeError('Run without -O.')
P = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('simplex_parent92',P.parent/'mw92_sections_20261115/verify.py')
parent92 = importlib.util.module_from_spec(spec); spec.loader.exec_module(parent92)
PIN = '52d3fea92df8eaa1c98776ce33e7a5fa42f32ca0'


def audit(export=None):
    if export:
        export.mkdir(parents=True,exist_ok=False)
    checkpoint = json.loads((P/'checkpoint.json').read_text())
    assert checkpoint['commit'] == PIN
    old = {a['dimension']:a for a in checkpoint['rows']}
    G, report = parent92.parent(); inv = G.inv()
    comparisons = parent92.rogers_comparison_intervals()
    rows = []
    for a in json.loads((P/'sections.json').read_text())['rows']:
        n = a['dimension']; assert n in (86,87,88,89)
        k = 92-n; B = fmpz_mat(a['normal_labels']); S = fmpz_mat(a['cluster_labels'])
        assert B.nrows() == k and B.ncols() == 92
        assert S.tolist() == [[0]*k]+[[int(i==j)for j in range(k)]for i in range(k)]
        A,T = B.transpose().hnf(transform=True)
        assert A == T*B.transpose() and abs(T.det()) == 1
        assert A.tolist()[:k] == [[int(i==j)for j in range(k)]for i in range(k)]
        lifts,N = fmpz_mat(T.tolist()[:k]),fmpz_mat(T.tolist()[k:])
        assert N*B.transpose() == fmpz_mat(n,k)
        assert lifts*B.transpose() == fmpz_mat(np.eye(k,dtype=int).tolist())
        K = fmpq_mat(B)*inv*fmpq_mat(B.transpose()); Q = N*G*N.transpose()
        assert Q.det() == G.det()*K.det()
        # These normals form a short, correlated simplex basis in the dual.
        assert all(str(K[i,i]) == '244/47' for i in range(k))
        shifts = fmpq_mat(S*lifts)-fmpq_mat(S)*K.inv()*fmpq_mat(B)*inv
        assert shifts*fmpq_mat(B.transpose()) == fmpq_mat(k+1,k)
        coords = shifts*fmpq_mat(G*N.transpose())*Q.inv()
        assert coords*fmpq_mat(N) == shifts
        beta = F(0)
        for i in range(k+1):
            for j in range(i):
                v = fmpq_mat([[int(S[i,t]-S[j,t])for t in range(k)]])
                beta = max(beta,F(str((v*K.inv()*v.transpose())[0,0])))
        assert 0 < beta < 16
        density = (k+1)**2*((16-beta)/4)**n/F(int(Q.det()))
        previous = F(old[n]['best_checked']['density_squared'])
        assert density > previous
        if n == 89:
            comparison = comparisons[n]
            lower,upper = map(F,(comparison['density_squared_lower'],comparison['density_squared_upper']))
            name = 'MW92 section using the full Rogers simplex theorem (derived interval)'
        else:
            lower = upper = {86:F(2**84,3),87:F(2**85),88:F(2**88)}[n]
            comparison = {'density_squared_lower':str(lower),'density_squared_upper':str(upper)}
            name = {86:'BCS88 section perpendicular to sqrt(2) A2',
                    87:'BCS88 section perpendicular to sqrt(2) A1',
                    88:'BCS88 E8/code construction'}[n]
        assert lower == F(old[n]['comparison']['density_squared'])
        assert density > upper
        rows.append({'dimension':n,'name':f'MW92 simplex antipode, {k+1} translates',
                     'kind':'new_explicit_periodic_packing','density_squared':str(density),
                     'center_density':parent92.decsqrt(density),'parent_dimension':92,
                     'normal_rank':k,'normal_gram':[[str(x)for x in r]for r in K.tolist()],
                     'normal_determinant':str(K.det()),'period_determinant':str(Q.det()),
                     'cosets':k+1,'cluster_squared_diameter':str(beta),
                     'minimum_squared_guarantee':str(16-beta),
                     'incremental_gain':parent92.gain_interval(density/previous),
                     'external_comparison_name':name,'external_comparison':comparison,
                     'external_gain':parent92.gain_interval(density/upper,density/lower),
                     'complete_integer_kernel_checked':True,
                     'all_projected_shifts_checked_exactly':True})
        if export:
            (export/f'dim{n}.json').write_text(json.dumps({
                'period_basis_in_parent':N.tolist(),'period_gram':Q.tolist(),
                'shifts_in_period_coordinates':[[str(x)for x in r]for r in coords.tolist()],
                'shifts_in_parent_coordinates':[[str(x)for x in r]for r in shifts.tolist()]},
                indent=2,default=int)+'\n')
    assert [a['dimension']for a in rows] == [86,87,88,89]
    return {'status':'verified_four_stronger_MW92_simplex_antipodes',
            'checkpoint_commit':PIN,'parent_audit':report,'rows':rows,
            'maximum_parent_or_target_dimension':92,
            'all_source_and_period_checks_exact':True,
            'minimum_attainment_claimed':False}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--write-result',action='store_true')
    ap.add_argument('--export',type=Path); args = ap.parse_args(); result = audit(args.export)
    if args.write_result: (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    else: assert json.loads((P/'results.json').read_text()) == result
    for a in result['rows']:
        print(a['dimension'],a['center_density'],'previous',a['incremental_gain']['display'],
              'external',a['external_gain']['display'])
