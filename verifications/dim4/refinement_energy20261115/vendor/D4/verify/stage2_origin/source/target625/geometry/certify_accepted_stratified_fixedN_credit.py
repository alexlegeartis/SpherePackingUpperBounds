"""Exact accepted-g surplus conditional on cardinality and a smaller radius.

The same unchanged accepted weighted angular polynomial is certified on a
smaller spherical domain. The radial hypothesis is explicit: all N neighbors
have norms in [1,R], so their mutual direction cosines are at most
1-1/(2R^2), which the exact proof cutoff must dominate.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import certify_fixedN_common_spherical as core
from certify_accepted_combined_fixedN_credit import common_record

RADII=[Q(x) for x in ['501473/500000','1005643/1000000','50379/50000','50767/50000','102597/100000']]


def verify(path,node_limit):
    data=json.loads(Path(path).read_text())
    common,exact=common_record(ROOT/data['accepted_combination']['accepted_origin_file'])
    assert data['accepted_combination']==exact
    assert list(map(Q,data['common_scalar_coefficients']))==list(map(Q,common['scalar_coefficients']))
    R=Q(data['condition_radial_radius']);s=Q(data['s'])
    assert R in RADII and 1<=R and R*R<2
    rho=1-Q(1)/(2*R*R)
    assert rho<=s<=Q(21,40)
    assert data['fixed_count'] in (23,24,25)
    report=core.verify(path,node_limit)
    scale=Q(exact['positive_normalization']);weightedB=scale*Q(report['B'])
    N=report['fixed_count'];credit=N*(Q(exact['weighted_accepted_B'])-weightedB)
    assert credit>0
    report.update(status='EXACT_RADIUS_STRATIFIED_FIXED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED',
        accepted_combination=exact,condition_radial_radius=str(R),condition_radius_index=RADII.index(R),
        geometric_cosine_upper=str(rho),proof_cosine_cutoff=str(s),
        geometric_cosine_proof='For1<=r,s<=R<sqrt(2), (r*r+s*s-1)/(2*r*s) is increasing in each radius; '
            'unit separation therefore gives direction cosine<=1-1/(2*R*R).',
        weighted_fixed_N_B=str(weightedB),weighted_fixed_N_B_decimal=float(weightedB),
        weighted_surplus_at_N=str(credit),weighted_surplus_decimal=float(credit),
        theorem='The unchanged accepted spherical sum has this surplus when its whole supported neighborhood '
            'has the stated exact countN and every point has norm<=the stated condition radius. '
            'This is a conditional theorem, not a claim that arbitrary neighborhoods satisfy that radius condition.',
        stratified_verifier_sha256=core.sha(Path(__file__)),
        accepted_combination_helper_sha256=core.sha(HERE/'certify_accepted_combined_fixedN_credit.py'),
        global_bound_proved=False,all_cardinality_claim=False)
    return report


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--epsilon',type=Q,required=True)
    p.add_argument('--pair-shift',type=Q,default=Q(1,1000000))
    p.add_argument('--radius',type=Q,required=True);p.add_argument('--cutoff',type=Q,required=True)
    p.add_argument('--node-limit',type=int,default=1200000);args=p.parse_args()
    folder=args.output.resolve();folder.mkdir(parents=True,exist_ok=False)
    common,exact=common_record(ROOT/'current/candidates/round9_six_graphs_degree18/origin.json')
    with np.load(args.source,allow_pickle=False) as source:
        assert float(source['s'])==float(args.cutoff)
        assert np.max(abs(source['scalar']-np.array([float(Q(x)) for x in common['scalar_coefficients']])))<1e-14
    common=dict(common,cosine_cutoff=str(args.cutoff))
    data=core.rationalize(args.source,common,args.epsilon,args.pair_shift)
    repair_file=core.BASE/'repair_s0525.json'
    data['repair']['gegenbauer_coefficients']=json.loads(repair_file.read_text())['gegenbauer_coefficients']
    data['accepted_combination']=exact;data['condition_radial_radius']=str(args.radius)
    data['numerical_source_sha256']=core.sha(args.source);data['repair_input_sha256']=core.sha(repair_file)
    path=folder/'certificate.json';path.write_text(json.dumps(data,indent=2)+'\n')
    try:report=verify(path,args.node_limit)
    except Exception as error:
        (folder/'NOT_ACCEPTED.json').write_text(json.dumps(dict(status='STRATIFIED_FIXED_N_SURPLUS_NOT_ACCEPTED',
            error=repr(error),certificate_sha256=core.sha(path)),indent=2)+'\n');raise
    (folder/'VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
