#!/usr/bin/env python3
"""Portable complete count audit. Preserve the three sibling source directories.
Rebuilds proof artifacts from frozen inputs, checks two exact cyclic replays and
an independent Sturm replay, tests -O rejection, and refreshes result/manifests.
No discovery scripts, external repository, network, or fixed workspace is used.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import hashlib
import json
import os
from pathlib import Path
from fractions import Fraction
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
PAIR=BASE/'density_kissing16_independent_pair_20261111'
FROZEN=BASE/'density_kissing16_reweight_frozen_20261111'
EXPECTED_INPUTS={'candidate.json':'b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573','factors.json':'6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954'}
EXPECTED_PROOF='d1d5333e9768677e53f924a0897358567e4164b7a09a6361046072572a669ca5'
EXPECTED_POLYNOMIAL='24085cd076fe5a8d9a4a5bda3b51a199e04b7ea9e82d7e922ce91a0d138bd3e7'

def require(condition,message):
    if not condition:raise RuntimeError(message)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def rel(path):return str(path.relative_to(BASE))

def main():
    start=time.monotonic();env=os.environ.copy()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):env[name]='1'
    env.pop('PYTHONOPTIMIZE',None)
    logs=HERE/'audit_logs';logs.mkdir(exist_ok=True)
    for name,digest in EXPECTED_INPUTS.items():require(sha(FROZEN/name)==digest,'Frozen input hash mismatch: '+name)
    # These are published mathematical proof artifacts and must not change.
    require(sha(HERE/'proof_events.json')==EXPECTED_PROOF,'Frozen traversal proof hash mismatch')
    require(sha(HERE/'polynomial.json')==EXPECTED_POLYNOMIAL,'Frozen polynomial hash mismatch')
    stages=[('pair_reconstruction',PAIR/'verify.py',[]),('pair_integer_replay',PAIR/'replay_chain.py',[]),('cyclic_construction',HERE/'bernstein_gram.py',['200000']),('cyclic_independent_construction_audit',PAIR/'audit_cyclic.py',[]),('cyclic_exact_integer_replay',HERE/'replay_exact.py',[]),('cyclic_independent_integer_replay',PAIR/'replay_gram_integer.py',[])]
    normal=[]
    for name,path,args in stages:
        t=time.monotonic();result=subprocess.run([sys.executable,str(path),*args],cwd=BASE,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (logs/(name+'.log')).write_text(result.stdout)
        require(result.returncode==0,f'{name} failed; inspect {rel(logs/(name+".log"))}')
        normal.append({'stage':name,'script':rel(path),'exit_code':result.returncode,'seconds':time.monotonic()-t})
        print('PASS',name,flush=True)
    for name,digest in EXPECTED_INPUTS.items():require(sha(FROZEN/name)==digest,'Frozen input changed: '+name)
    require(sha(HERE/'proof_events.json')==EXPECTED_PROOF,'Mathematical proof events changed')
    require(sha(HERE/'polynomial.json')==EXPECTED_POLYNOMIAL,'Mathematical polynomial changed')
    # All public Python files in the two audit directories, including imported
    # proof entrypoints, reject optimized Python before executing any checks.
    rejected=[]
    for path in sorted(PAIR.glob('*.py'))+sorted(HERE.glob('*.py')):
        result=subprocess.run([sys.executable,'-O',str(path)],cwd=BASE,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        require(result.returncode!=0 and 'Certificate verification requires assertions' in result.stdout,'-O did not trigger explicit guard: '+rel(path))
        rejected.append({'script':rel(path),'exit_code':result.returncode,'explicit_optimization_guard':True})
    # -OO is also checked for the top-level entrypoint.
    result=subprocess.run([sys.executable,'-OO',str(HERE/'audit_count.py')],cwd=BASE,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    require(result.returncode!=0 and 'Certificate verification requires assertions' in result.stdout,'-OO did not trigger explicit guard')
    print('PASS optimized-mode rejection for',len(rejected),'entrypoints and top-level -OO',flush=True)
    pair=read(PAIR/'result.json');sturm=read(PAIR/'replay_result.json');audit=read(PAIR/'cyclic_audit_result.json');exact=read(HERE/'exact_replay_result.json');independent=read(PAIR/'gram_integer_replay_result.json');cover=read(HERE/'bernstein_gram_result.json')
    require(pair['all_source_lambdas_positive'] and pair['all_theta_a_r_nonnegative'],'PSD/coefficient gate missing')
    require(pair['dimension']==16 and pair['cap_exact']=='5003/10000','Wrong count problem')
    require(sturm['status']=='PASS' and all(x['roots']==0 and x['all_remainder_identities_verified'] for x in sturm['certificates']),'Sturm replay gate missing')
    require(audit['independent_factor_form_matches_exact_cyclic_polynomial'] and audit['exact_permutation_symmetry'],'Independent polynomial gate missing')
    require(exact['verified'] and exact['covered_unit_cube_volume_exact']=='1','Exact cyclic replay incomplete')
    require(independent['status'].startswith('PASS') and independent['full_unit_cube_covered_exact']=='1','Independent integer replay incomplete')
    require(cover['verified'] and exact['events']==independent['checked']==cover['checked']==2155,'Cover/event inconsistency')
    U=Fraction(pair['objective_exact']);count=U.numerator//U.denominator;require(count==7802,'Unexpected final count')
    certificate={'status':'PASS: complete positive-width spherical counting witness','dimension':16,'cap_exact':'5003/10000','spherical_count_integer_upper_bound':count,'real_upper_bound_exact':str(U),'real_upper_bound_enclosure':pair['objective_enclosure'],'input_sha256':EXPECTED_INPUTS,'psd_and_nonnegative_coefficients_pass':True,'pair_uniform_upper_bound_exact':'-9/10000000','pair_proof':'Exact rational Sturm chains, independently replayed with Python integers/Fraction','cyclic_uniform_upper_bound_exact':'-'+independent['rigorous_uniform_cyclic_negative_margin_exact'],'cyclic_proof':'Complete Bernstein covering; independent exact polynomial reconstruction; two full integer acceptance replays','covered_normalized_cube_volume_exact':'1','cyclic_checked':exact['events'],'cyclic_counts':exact['counts'],'source_contact_SOS_replay_needed':False,'source_contact_SOS_replay_performed':False,'density_claim':False,'density_scope':'Only the spherical counting gate is certified here. Density radial/profile/geometric gates require separate verification.','portable_top_level_audit':'audit_count.py','source_directory_dependencies':[FROZEN.name,PAIR.name,HERE.name],'artifacts_sha256':{rel(p):sha(p) for p in [FROZEN/'candidate.json',FROZEN/'factors.json',PAIR/'verify.py',PAIR/'replay_chain.py',PAIR/'audit_cyclic.py',PAIR/'replay_gram_integer.py',PAIR/'result.json',PAIR/'replay_result.json',PAIR/'cyclic_audit_result.json',PAIR/'gram_integer_replay_result.json',PAIR/'gram_replayed_leaves.json',HERE/'bernstein.py',HERE/'bernstein_gram.py',HERE/'polynomial.json',HERE/'proof_events.json',HERE/'replay_exact.py',HERE/'exact_replay_result.json',HERE/'audit_count.py']}}
    (HERE/'CERTIFIED_COUNT.json').write_text(json.dumps(certificate,indent=2)+'\n')
    output={'status':'PASS','purpose':'Portable full spherical-count audit, not log-only validation','spherical_count_integer_upper_bound':count,'normal_stages':normal,'optimization_rejections':rejected,'top_level_double_optimization_rejected':True,'frozen_inputs_unchanged':True,'mathematical_proof_events_unchanged':True,'mathematical_polynomial_unchanged':True,'source_directory_dependencies':[FROZEN.name,PAIR.name,HERE.name],'package_dependencies':['Python >=3.11','numpy','python-flint'],'tested_python':sys.version,'seconds':time.monotonic()-start}
    (HERE/'AUDIT_RESULT.json').write_text(json.dumps(output,indent=2)+'\n')
    # No manifest includes itself, another manifest, bytecode, or subdirectories.
    for directory in (PAIR,HERE):
        files=[p for p in directory.iterdir() if p.is_file() and p.suffix != '.log' and p.name not in ('SHA256SUMS','PACKAGE_SHA256SUMS')]
        (directory/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(files)))
    files=[p for directory in (FROZEN,PAIR,HERE) for p in directory.iterdir() if p.is_file() and p.suffix != '.log' and p.name not in ('SHA256SUMS','PACKAGE_SHA256SUMS')]
    (HERE/'PACKAGE_SHA256SUMS').write_text(''.join(f'{sha(p)}  {rel(p)}\n' for p in sorted(files)))
    print('PASS complete count certificate:',count,'points at cap 5003/10000',flush=True)

if __name__=='__main__':main()
