"""Fresh exact bound H_g(N)<=N*B_N for one stated N and common rational g.

The fixed-N anchor quadratic is certified by integer Gram rays after the
substitution radial_count=N-1. A symbolic PSD repair preserves the common g.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np
import fixed_cardinality_spherical_algebra as algebra

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE=ROOT/'optimization/continuation_20261018/base/new_work/spherical_sdp'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def rationalize(source, common, epsilon, pair_shift, feature_scale=10**9):
    with np.load(source,allow_pickle=False) as z:
        matrices=z['matrices'].copy();v=z['v'].copy();kappa=float(z['kappa'])
        N=int(z['fixed_count']);d=int(z['radial_degree']);L=int(z['angular_degree'])
        assert str(z['radial_basis'])=='chebyshev' and N in (23,24,25)
    g=list(map(Q,common['scalar_coefficients']));n=N-1
    reduced=np.block([[np.array([[kappa+2*n*v[0]+n*n*matrices[0][0,0]]]),
                       (v[1:]+n*matrices[0][0,1:])[None,:]],
                      [(v[1:]+n*matrices[0][0,1:])[:,None],matrices[0][1:,1:]]])
    blocks=[]
    for matrix in [reduced]+list(matrices[1:]):
        eigenvalues,eigenvectors=np.linalg.eigh((matrix+matrix.T)/2)
        factors=np.rint(eigenvectors*np.sqrt(np.maximum(eigenvalues,0))[None,:]*feature_scale)
        blocks.append([[int(factors[i,j]) for i in range(len(matrix))] for j in range(len(matrix))])
    scale=feature_scale**2
    free=[int(round(float(x)*scale)) for x in [v[0],matrices[0][0,0],*matrices[0][0,1:]]]
    denominator=math.lcm(*(a.denominator for a in g))
    repair=json.loads((BASE/'repair_s0515.json').read_text())
    return dict(dimension=4,s=common['cosine_cutoff'],radial_degree=d,radial_basis='chebyshev',
        angular_degree=L,scalar_degree=len(g)-1,feature_scale=feature_scale,
        scalar_scale=denominator,scalar=[int(a*denominator) for a in g],blocks=blocks,
        fixed_count=N,fixed_origin_free=free,verification_goal='fixed_cardinality_signed_energy',
        common_scalar_coefficients=common['scalar_coefficients'],
        repair=dict(scale=str(epsilon),pair_shift=str(pair_shift),
                    gegenbauer_coefficients=repair['gegenbauer_coefficients']),
        source_status='NUMERICAL_FIXED_N_PROPOSAL_REQUIRING_FRESH_EXACT_PROOF')

def verify(path,node_limit):
    started=time.time();data=json.loads(Path(path).read_text());N=data['fixed_count'];n=N-1
    triangle,td,pair,pd,gram,B,s=algebra.polynomials(data)
    algebra.audit_identities(data,triangle,td,pair,pd)
    g=[Q(a,data['scalar_scale']) for a in data['scalar']]
    assert list(map(str,g))==list(map(str,map(Q,data['common_scalar_coefficients'])))
    repair=data['repair'];epsilon=Q(repair['scale']);eta=Q(repair['pair_shift'])
    coefficients=list(map(Q,repair['gegenbauer_coefficients']))
    assert epsilon>0 and eta>=0 and coefficients[0]>0 and all(a>=0 for a in coefficients)
    polynomial=[Q(0)]*len(coefficients)
    for ell,a in enumerate(coefficients):
        for j in range(ell//2+1):
            polynomial[ell-2*j]+=a*Q((-1)**j*math.comb(ell-j,j)*2**(ell-2*j),ell+1)
    rd=math.lcm(*(a.denominator for a in polynomial))
    rc=np.array([int(a*rd) for a in polynomial],dtype=object)
    repair_sign=algebra.certify(rc,rd,[(Q(-1),s)],node_limit=node_limit)
    assert repair_sign['status']=='CERTIFIED',repair_sign
    H=sum(coefficients)/coefficients[0]-1;assert H>0
    delta=epsilon*H+eta
    # S'_anchor=S_anchor+epsilon*sum_ij T(t_ij)
    #             -delta*(sum_i1-(N-1)).
    # The last term is identically0 at exact N; T is PSD with diagonalH and
    # distinct entries<=-1. Thus common g is unchanged, pair decreaseseta,
    # triangle decreasesat least3epsilon, and kappa growsby(N-1)delta.
    new_td=math.lcm(td,(3*epsilon).denominator)
    triangle*=new_td//td;triangle[0,0,0]-=int(3*epsilon*new_td);td=new_td
    new_pd=math.lcm(pd,eta.denominator)
    pair*=new_pd//pd;pair[0]-=int(eta*new_pd);pd=new_pd
    pair_result=algebra.certify(pair,pd,[(Q(-1),s)],node_limit=node_limit)
    print('PAIR',N,pair_result,flush=True);assert pair_result['status']=='CERTIFIED'
    triangle_result=algebra.certify(triangle,td,[(Q(-1),s)]*3,gram=gram,node_limit=node_limit)
    print('TRIANGLE',N,triangle_result,flush=True);assert triangle_result['status']=='CERTIFIED'
    finalB=B+n*delta
    report=dict(status='EXACT_FIXED_CARDINALITY_SIGNED_SPHERICAL_ENERGY_VERIFIED',
        fixed_count=N,cosine_cutoff=str(s),common_scalar_coefficients=list(map(str,g)),
        B=str(finalB),B_decimal=float(finalB),HN_upper=str(N*finalB),
        raw_B=str(B),repair=dict(scale=str(epsilon),pair_shift=str(eta),diagonal=str(H),
            count_identity_delta=str(delta),B_increment=str(n*delta),polynomial_sign=repair_sign),
        pair=pair_result,triangle=triangle_result,
        exact_positivity='Integer Gram rays for the N-1 radial-count restriction plus PSD T; common scalar g unchanged.',
        all_cardinality_claim=False,global_bound_proved=False,
        certificate_sha256=sha(path),verifier_sha256=sha(Path(__file__)),
        algebra_sha256=sha(HERE/'fixed_cardinality_spherical_algebra.py'),seconds=time.time()-started)
    return report

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--epsilon',type=Q,required=True)
    p.add_argument('--pair-shift',type=Q,default=Q(1,10000000))
    p.add_argument('--node-limit',type=int,default=1200000)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    common_file=HERE/'fresh_s0515_degree12_exact_eps3_extended/verified_local_energy_ray.json'
    common=json.loads(common_file.read_text())['records'][0]
    data=rationalize(args.source,common,args.epsilon,args.pair_shift)
    data['common_source_sha256']=sha(common_file);data['numerical_source_sha256']=sha(args.source)
    certificate=args.output/'certificate.json';certificate.write_text(json.dumps(data,indent=2)+'\n')
    try:report=verify(certificate,args.node_limit)
    except Exception as error:
        (args.output/'NOT_ACCEPTED.json').write_text(json.dumps(dict(status='EXACT_FIXED_N_NOT_ACCEPTED',error=repr(error),certificate_sha256=sha(certificate)),indent=2)+'\n')
        raise
    (args.output/'VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
