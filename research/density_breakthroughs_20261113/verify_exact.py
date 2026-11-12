#!/usr/bin/env python3
"""Replay the compact exact claims; keep numerical screens out of this gate."""
from pathlib import Path
import json,os,subprocess,sys,time
P=Path(__file__).resolve().parent
scripts=[
    'contact_walk/verify36_local.py',
    'mixed_codes/verify.py',
    'sections/audit_n17.py',
    'catalogue_screen/verify56_pencil.py',
    'gluing/verify_extension_thresholds.py',
    'gluing/verify_normal_witnesses.py',
    'sections/verify_finite_results.py',
    'seed_deformations/audit_seed_results.py',
]
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
rows=[]
for script in scripts:
    start=time.time()
    run=subprocess.run([sys.executable,str(P/script)],cwd=P.parents[1],env=env,capture_output=True,text=True)
    row={'script':script,'passed':run.returncode==0,'seconds':time.time()-start}
    rows.append(row);print(json.dumps(row),flush=True)
    if run.returncode:
        print(run.stdout);print(run.stderr,file=sys.stderr)
        (P/'exact_verification_summary.json').write_text(json.dumps({'checks':rows,'all_passed':False},indent=2))
        raise SystemExit(run.returncode)
report={'all_passed':True,'checks':rows,'new_density_record':False,'scope':'Exact finite certificates and comparisons. Prior published ambient minima and the existing36D minimum certificate remain explicit inputs; numerical contact-cone screens are excluded.'}
(P/'exact_verification_summary.json').write_text(json.dumps(report,indent=2))
