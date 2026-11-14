"""Materialize the exact successful radial witnesses with complete replay tools."""
from pathlib import Path
import json,shutil,hashlib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
old={45:'dim45',46:'dim46/refinement_extended_pair',47:'dim47/refinement_degree47'}
degrees={45:47,46:63,47:47}
strict={45:'0.000001287489137',46:'0.0000008756454737',47:'0.0000005972403006'}
previous={45:'0.000001311423575',46:'0.0000008777229084',47:'0.0000005998306274'}
angular_audit=(ROOT/'verifications/dim46/refinement_extended_pair/independent/audit_angular_cyclic.py').read_text()
angular_audit=angular_audit.replace("D/'inputs/angular46.json'","D/'inputs/angular.json'")
angular_audit=angular_audit.replace("assert n == 46 and r['degree'] == 19 and r['alpha'] == 40",
    "cfg=json.loads((D/'inputs/config.json').read_text())\nassert n == cfg['dimension'] and n in [45,46,47] and r['degree'] == 19 and r['alpha'] == 40")
for n,degree in degrees.items():
    package=f'verifications/dim{n}/refinement_radial20261114'
    D=ROOT/package;D.mkdir(parents=True,exist_ok=True)
    S=D/'continuum_three_point';S.mkdir(exist_ok=True)
    I=D/'independent';I.mkdir(exist_ok=True)
    (I/'inputs').mkdir(exist_ok=True)
    src=ROOT/'verifications'/old[n]/'continuum_three_point'
    for file in ['cyclic_angular40_exact.json','independent_cyclic40_angular.py','independent_cyclic40_domain.py']:
        shutil.copy2(src/file,S/file)
    candidate=HERE/f'n{n}_d{degree}_c40'
    data=json.loads((candidate/'radial_certificate.json').read_text())
    assert data['radial_checks_pass'] and data['M_guard_exact']=='1000001/1000000'
    keys=['dimension','degree','physical_polynomial','fourier_polynomial','M_exact','B_exact',
          'b0_exact','M_guard_exact','objective_exact','coupled_upper_X_exact','physical_tail_X_exact']
    (D/'certificate.json').write_text(json.dumps({k:data[k] for k in keys},indent=2)+'\n')
    cfg={'dimension':n,'package':package,'strict_density_upper_bound':strict[n],
         'previous_repository_strict_upper':previous[n],'previous_repository_package':'verifications/'+old[n],
         'radial_correction_degree':degree,'minimum_squared_distance':'6','kappa_exact':'17/8',
         'base_commit':'72028bb0c5eaaaeede50cd418ca09339ee404285',
         'angular_polynomial_unchanged':True,'external_neighbor_count_inputs':False}
    (D/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    for file in ['certificate.json','config.json']:
        shutil.copy2(D/file,I/'inputs'/file)
    shutil.copy2(S/'cyclic_angular40_exact.json',I/'inputs/angular.json')
    shutil.copy2(HERE/'verify_radial_primary.py',D/'verify_radial_primary.py')
    shutil.copy2(HERE/'audit_radial_independent.py',I/'audit_radial.py')
    (I/'audit_angular_cyclic.py').write_text(angular_audit)
    (D/'requirements.txt').write_text('python-flint==0.9.0\n')
    runner='''"""Replay every primary and independent gate of this frozen certificate."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
D=Path(__file__).resolve().parent
cfg=json.loads((D/'config.json').read_text())
for name,want in json.loads((D/'source_hashes.json').read_text()).items():
    assert hashlib.sha256((D/name).read_bytes()).hexdigest()==want,name
rows=json.loads((D.parents[2]/'RESULTS.json').read_text())
row=next(r for r in rows if r['dimension']==cfg['dimension'])
assert row['package']==cfg['package']
assert Fraction(row['strict_density_upper_bound'])==Fraction(cfg['strict_density_upper_bound'])
logs=[]
for script in ['continuum_three_point/independent_cyclic40_angular.py',
               'continuum_three_point/independent_cyclic40_domain.py',
               'verify_radial_primary.py','independent/audit_angular_cyclic.py',
               'independent/audit_radial.py']:
    result=subprocess.run([sys.executable,str(D/script)],cwd=D,capture_output=True,text=True)
    log=script.replace('/','_').removesuffix('.py')+'.log'
    (D/log).write_text(result.stdout+result.stderr)
    print(result.stdout,flush=True)
    assert result.returncode==0,(script,result.stderr)
    logs.append({'script':script,'returncode':result.returncode})
out={'pass':True,'dimension':cfg['dimension'],'strict_density_upper_bound':cfg['strict_density_upper_bound'],
     'complete_primary_and_independent_replay':True,'executed':logs,
     'angular_data_unchanged':True,'arithmetic_backend':'FLINT/Arb shared; independent proof implementations'}
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2)+'\\n')
print('FULL CERTIFICATE PASS',cfg['dimension'],cfg['strict_density_upper_bound'],flush=True)
'''
    (D/'verify.py').write_text(runner)
    files=['certificate.json','config.json','verify_radial_primary.py',
           'continuum_three_point/cyclic_angular40_exact.json',
           'continuum_three_point/independent_cyclic40_angular.py',
           'continuum_three_point/independent_cyclic40_domain.py',
           'independent/inputs/certificate.json','independent/inputs/config.json',
           'independent/inputs/angular.json','independent/audit_radial.py','independent/audit_angular_cyclic.py']
    manifest={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in files}
    (D/'source_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(package,strict[n])
