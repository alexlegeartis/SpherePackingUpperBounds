"""Replay every primary and independent gate of this frozen certificate."""
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
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('FULL CERTIFICATE PASS',cfg['dimension'],cfg['strict_density_upper_bound'],flush=True)
