"""Aggregate the four separately audited improvements and check their checkpoint."""
import argparse,json
from fractions import Fraction as F
from common import P,checkpoint
from verify_p48 import setup,audit as audit_p48,D
from verify_section53 import audit as audit53
from verify_adaptive import audit as audit_adaptive
from verify_taylor_covering import audit as audit_covering
if not __debug__:raise RuntimeError('Run without -O.')

def audit():
    files,rows,Q=setup()
    reports={'p48_results.json':audit_p48(D/'replay',files,rows,Q),
             'section53_results.json':audit53(),'adaptive_results.json':audit_adaptive(),
             'covering_results.json':audit_covering()}
    rerun=audit_p48(D/'replay_rerun',files,rows,Q)
    assert rerun==reports['p48_results.json']
    out=[]
    for name,r in reports.items():
        assert json.loads((P/name).read_text())==r,(name,'Stored report stale')
        for a in r['rows']:
            out.append({**a,'proof_report':name})
    out.sort(key=lambda a:a['dimension'])
    assert [a['dimension']for a in out]==[50,51,52,53,57,77,81,82,85,93,94,95,96]
    for a in out:assert F(a['density_squared'])>F(checkpoint()[a['dimension']]['density_squared'])
    return {'status':'13_strictly_stronger_certified_lower_bounds','maximum_parent_or_target_dimension':96,
            'checkpoint_commit':json.loads((P/'checkpoint.json').read_text())['commit'],'rows':out,
            'all_four_exact_audits_passed':True,'fresh_full_P48_rerun_matches_primary':True,
            'fresh_full_P48_rerun_completed_nodes':rerun['completed_nodes'],
            'proof_dependencies':'Published parent minima and computational Grams, standard modular theta and binary quadratic-space facts, Coxeter-Few-Rogers theorem; P48 enumeration arithmetic assumptions stated explicitly.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args();r=audit()
    if args.write_result:(P/'results.json').write_text(json.dumps(r,indent=2,default=int)+'\n')
    else:assert json.loads((P/'results.json').read_text())==r
    for a in r['rows']:print(a['dimension'],a['center_density'],a['gain_over_checkpoint_percent'])
