"""Fresh exact fixed-N surplus in the unchanged accepted spherical sum.

The numerical common g is rescaled to exactly equal the accepted rational
weighted sum. Integer Gram rays and the common-g preserving harmonic repair
are delegated to the already independently audited fixed-N checker.
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


def common_record(origin):
    records=json.loads(Path(origin).read_text())['appended_spherical_local_energy_rays']
    assert len(records)==2
    assert all(r['status']=='EXACT_SPHERICAL_LOCAL_ENERGY_RAY' and
               Q(r['cosine_cutoff'])==Q(21,40) and Q(r['R'])==Q(102597,100000)
               and Q(r['global_origin'])==0 for r in records)
    weights=[Q(r['weight']) for r in records]
    assert min(weights)>0
    scale=sum(weights)
    g=[sum(w*Q(r['scalar_coefficients'][j]) for w,r in zip(weights,records)) for j in range(13)]
    B=sum(w*Q(r['B']) for w,r in zip(weights,records))
    normalized=[v/scale for v in g]
    result=dict(cosine_cutoff='21/40',scalar_coefficients=list(map(str,normalized)))
    exact=dict(accepted_origin_file=str(Path(origin).relative_to(ROOT)),
        accepted_origin_sha256=core.sha(origin),weights=list(map(str,weights)),
        positive_normalization=str(scale),weighted_accepted_g=list(map(str,g)),
        weighted_accepted_B=str(B),normalized_accepted_B=str(B/scale),
        accepted_source_sha256={r['source_certificate']:r['source_sha256'] for r in records})
    return result,exact


def verify(certificate,node_limit):
    path=Path(certificate);data=json.loads(path.read_text())
    origin=ROOT/data['accepted_combination']['accepted_origin_file']
    common,exact=common_record(origin)
    assert data['accepted_combination']==exact
    assert list(map(Q,data['common_scalar_coefficients']))==list(map(Q,common['scalar_coefficients']))
    assert data['s']=='21/40' and data['fixed_count'] in (23,24,25)
    report=core.verify(path,node_limit)
    scale=Q(exact['positive_normalization']);oldB=Q(exact['weighted_accepted_B'])
    newB=scale*Q(report['B']);N=report['fixed_count']
    credit=N*(oldB-newB)
    assert credit>0,'No useful surplus in the unchanged accepted spherical energy.'
    report.update(status='EXACT_FIXED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED',
        accepted_combination=exact,normalized_fixed_N_B=report['B'],
        weighted_fixed_N_B=str(newB),weighted_fixed_N_B_decimal=float(newB),
        weighted_surplus_at_N=str(credit),weighted_surplus_decimal=float(credit),
        theorem='For every admissible spherical code of the stated exact cardinalityN, '
                'the unchanged accepted sum of two spherical local energies is at least '
                'the stated positive weighted surplus. No claim for other N is made here.',
        accepted_combination_verifier_sha256=core.sha(Path(__file__)),
        core_verifier_sha256=core.sha(HERE/'certify_fixedN_common_spherical.py'),
        global_bound_proved=False,all_cardinality_claim=False)
    return report


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--epsilon',type=Q,required=True)
    parser.add_argument('--pair-shift',type=Q,default=Q(1,1000000))
    parser.add_argument('--node-limit',type=int,default=1200000)
    args=parser.parse_args();folder=args.output.resolve();folder.mkdir(parents=True,exist_ok=False)
    origin=ROOT/'current/candidates/round9_six_graphs_degree18/origin.json'
    common,exact=common_record(origin)
    with np.load(args.source,allow_pickle=False) as data:
        assert float(data['s'])==.525
        assert np.max(abs(data['scalar']-np.array([float(Q(v)) for v in common['scalar_coefficients']])))<1e-14
    data=core.rationalize(args.source,common,args.epsilon,args.pair_shift)
    repair_file=core.BASE/'repair_s0525.json'
    repair=json.loads(repair_file.read_text());assert repair['s']=='21/40'
    data['repair']['gegenbauer_coefficients']=repair['gegenbauer_coefficients']
    data['accepted_combination']=exact
    data['numerical_source_file']=str(args.source.resolve().relative_to(ROOT))
    data['numerical_source_sha256']=core.sha(args.source)
    data['repair_input_sha256']=core.sha(repair_file)
    certificate=folder/'certificate.json'
    certificate.write_text(json.dumps(data,indent=2)+'\n')
    try:report=verify(certificate,args.node_limit)
    except Exception as error:
        (folder/'NOT_ACCEPTED.json').write_text(json.dumps(dict(status='EXACT_ACCEPTED_SPHERICAL_SURPLUS_NOT_ACCEPTED',
            error=repr(error),certificate_sha256=core.sha(certificate)),indent=2)+'\n')
        raise
    (folder/'VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
