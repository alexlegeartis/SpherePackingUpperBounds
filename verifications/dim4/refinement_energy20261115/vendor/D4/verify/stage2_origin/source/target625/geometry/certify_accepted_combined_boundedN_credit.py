"""Exact common-g surplus valid for every spherical count 0<=N<=M.

An unrestricted integer-Gram anchor S is repaired by
  epsilon sum_ij T(t_ij) + delta*((M-1)-n), n=N-1,
where T is positive definite, T(1)=H, T(t)<=-1 off diagonal,
and delta=epsilon*H+eta. The repair is nonnegative for n<=M-1,
keeps the common scalar g, and uses one count-independent budget.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import json
import math
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
BASE=ROOT/'optimization/continuation_20261018/base/new_work/spherical_sdp'
sys.path.insert(0,str(HERE))
from certify_accepted_combined_fixedN_credit import common_record
from certify_fixedN_common_spherical import sha
sys.path.insert(0,str(BASE))
import certify as algebra

ALGEBRA_SHA='aeb873edff74112cd0dd0293d1c77be91aec7ec405aa1013b9b92088db46fca5'


def rationalize(source,common,epsilon,eta,maximum_count):
    with np.load(source,allow_pickle=False) as z:
        assert int(z['fixed_count'])==0 and str(z['radial_basis'])=='chebyshev'
        assert float(z['s'])==.525
        matrices=z['matrices'].copy();v=z['v'].copy();kappa=float(z['kappa'])
        d=int(z['radial_degree']);L=int(z['angular_degree'])
        # Equality-constrained SDP variables have numerical residuals. The
        # exact scalar below is replaced by the unchanged accepted rational
        # g, and every resulting sign inequality is freshly replayed.
        scalar_deviation=float(np.max(abs(z['scalar']-np.array([float(Q(v)) for v in common['scalar_coefficients']]))))
        assert scalar_deviation<1e-6
    augmented=np.block([[np.array([[kappa]]),v[None,:]],[v[:,None],matrices[0]]])
    feature_scale=10**9;blocks=[]
    for matrix in [augmented]+list(matrices[1:]):
        eigenvalues,eigenvectors=np.linalg.eigh((matrix+matrix.T)/2)
        factors=np.rint(eigenvectors*np.sqrt(np.maximum(eigenvalues,0))[None,:]*feature_scale)
        blocks.append([[int(factors[i,j]) for i in range(len(matrix))] for j in range(len(matrix))])
    g=list(map(Q,common['scalar_coefficients']));den=math.lcm(*(x.denominator for x in g))
    repair=json.loads((BASE/'repair_s0525.json').read_text())
    return dict(dimension=4,s='21/40',radial_degree=d,radial_basis='chebyshev',
        angular_degree=L,scalar_degree=len(g)-1,feature_scale=feature_scale,
        scalar_scale=den,scalar=[int(x*den) for x in g],blocks=blocks,
        verification_goal='local_energy',common_scalar_coefficients=list(map(str,g)),
        numerical_scalar_equality_residual=scalar_deviation,
        bounded_cardinality_common_g_repair=dict(maximum_code_count=maximum_count,
            scale=str(epsilon),pair_shift=str(eta),gegenbauer_coefficients=repair['gegenbauer_coefficients']))


def verify(path,node_limit):
    started=time.time();assert sha(BASE/'certify.py')==ALGEBRA_SHA
    data=json.loads(Path(path).read_text());assert not data.get('fixed_count')
    assert data['verification_goal']=='local_energy' and data['s']=='21/40'
    common,exact=common_record(ROOT/data['accepted_combination']['accepted_origin_file'])
    assert data['accepted_combination']==exact
    assert list(map(Q,data['common_scalar_coefficients']))==list(map(Q,common['scalar_coefficients']))
    triangle,td,pair,pd,gram,B,s=algebra.polynomials(data)
    algebra.audit_identities(data,triangle,td,pair,pd)
    assert [Q(x,data['scalar_scale']) for x in data['scalar']]==list(map(Q,common['scalar_coefficients']))
    repair=data['bounded_cardinality_common_g_repair']
    M=repair['maximum_code_count'];assert type(M) is int and 1<=M<=1000
    epsilon=Q(repair['scale']);eta=Q(repair['pair_shift'])
    coefficients=list(map(Q,repair['gegenbauer_coefficients']))
    assert epsilon>0 and eta>=0 and coefficients[0]>0 and all(c>=0 for c in coefficients)
    polynomial=[Q(0)]*len(coefficients)
    for ell,c in enumerate(coefficients):
        for j in range(ell//2+1):
            polynomial[ell-2*j]+=c*Q((-1)**j*math.comb(ell-j,j)*2**(ell-2*j),ell+1)
    rd=math.lcm(*(c.denominator for c in polynomial))
    rc=np.array([int(c*rd) for c in polynomial],dtype=object)
    repair_sign=algebra.certify(rc,rd,[(Q(-1),s)],node_limit=node_limit)
    assert repair_sign['status']=='CERTIFIED',repair_sign
    H=sum(coefficients)/coefficients[0]-1;assert H>0
    delta=epsilon*H+eta
    new_td=math.lcm(td,(3*epsilon).denominator)
    triangle*=new_td//td;triangle[0,0,0]-=int(3*epsilon*new_td);td=new_td
    new_pd=math.lcm(pd,eta.denominator)
    pair*=new_pd//pd;pair[0]-=int(eta*new_pd);pd=new_pd
    pair_result=algebra.certify(pair,pd,[(Q(-1),s)],node_limit=node_limit)
    print('PAIR',pair_result,flush=True);assert pair_result['status']=='CERTIFIED'
    triangle_result=algebra.certify(triangle,td,[(Q(-1),s)]*3,gram=gram,node_limit=node_limit)
    print('TRIANGLE',triangle_result,flush=True);assert triangle_result['status']=='CERTIFIED'
    finalB=B+(M-1)*delta
    scale=Q(exact['positive_normalization']);weightedB=scale*finalB
    perN=Q(exact['weighted_accepted_B'])-weightedB
    assert perN>0
    return dict(status='EXACT_BOUNDED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED',
        maximum_code_count=M,cosine_cutoff=str(s),common_scalar_coefficients=common['scalar_coefficients'],
        raw_B=str(B),normalized_B=str(finalB),normalized_B_decimal=float(finalB),
        weighted_B=str(weightedB),weighted_B_decimal=float(weightedB),
        weighted_surplus_per_neighbor=str(perN),weighted_surplus_per_neighbor_decimal=float(perN),
        weighted_surplus_at_N24=str(24*perN),weighted_surplus_at_N24_decimal=float(24*perN),
        accepted_combination=exact,
        repair=dict(scale=str(epsilon),pair_shift=str(eta),diagonal=str(H),
            count_slack_coefficient=str(delta),B_increment=str((M-1)*delta),polynomial_sign=repair_sign),
        pair=pair_result,triangle=triangle_result,
        exact_positivity='Unrestricted integer-Gram anchor plus epsilon*sum(T) and delta*((M-1)-n); '
                        'both added terms nonnegative for all code counts1<=N<=M; N0 is trivial.',
        cardinality_scope=f'Every integer code size0<=N<={M}; no statement that all codes satisfy this count bound.',
        common_g_unchanged=True,global_bound_proved=False,
        certificate_sha256=sha(path),verifier_sha256=sha(Path(__file__)),
        algebra_sha256=ALGEBRA_SHA,
        accepted_combination_helper_sha256=sha(HERE/'certify_accepted_combined_fixedN_credit.py'),
        seconds=time.time()-started)


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--epsilon',type=Q,required=True)
    p.add_argument('--pair-shift',type=Q,default=Q(1,1000000));p.add_argument('--maximum-count',type=int,default=26)
    p.add_argument('--node-limit',type=int,default=1200000);args=p.parse_args()
    folder=args.output.resolve();folder.mkdir(parents=True,exist_ok=False)
    common,exact=common_record(ROOT/'current/candidates/round9_six_graphs_degree18/origin.json')
    data=rationalize(args.source,common,args.epsilon,args.pair_shift,args.maximum_count)
    data['accepted_combination']=exact;data['numerical_source_sha256']=sha(args.source)
    path=folder/'certificate.json';path.write_text(json.dumps(data,indent=2)+'\n')
    try:report=verify(path,args.node_limit)
    except Exception as error:
        (folder/'NOT_ACCEPTED.json').write_text(json.dumps(dict(status='BOUNDED_N_SURPLUS_NOT_ACCEPTED',
            error=repr(error),certificate_sha256=sha(path)),indent=2)+'\n');raise
    (folder/'VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
