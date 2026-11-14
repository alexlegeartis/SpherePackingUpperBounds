"""Replay every gate of this frozen certificate; assertions are required."""
import hashlib,json,subprocess,sys
from fractions import Fraction
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions must remain enabled')
D=Path(__file__).resolve().parent
cfg=json.loads((D/'config.json').read_text())
for name,want in json.loads((D/'source_hashes.json').read_text()).items():
    assert hashlib.sha256((D/name).read_bytes()).hexdigest()==want,name
row=next(r for r in json.loads((D.parents[2]/'RESULTS.json').read_text()) if r['dimension']==cfg['dimension'])
assert Fraction(row['strict_density_upper_bound'])<=Fraction(cfg['strict_density_upper_bound'])
assert row['package']==cfg['package'] or cfg['package'] in row.get('earlier_repository_packages',[])
executed=[]
for script in ['primary_angular.py','primary_radial.py','independent_angular.py','independent_radial.py']:
    result=subprocess.run([sys.executable,str(D.parents[1]/'profile_enrichment_20261114'/script),str(D)],capture_output=True,text=True)
    (D/(script.removesuffix('.py')+'.log')).write_text(result.stdout+result.stderr)
    print(result.stdout,flush=True)
    assert result.returncode==0,(script,result.stderr)
    executed.append({'script':script,'returncode':result.returncode})
out={'pass':True,'dimension':cfg['dimension'],'strict_density_upper_bound':cfg['strict_density_upper_bound'],
     'complete_primary_and_independent_replay':True,'executed':executed,
     'source_manifest_sha256':hashlib.sha256((D/'source_hashes.json').read_bytes()).hexdigest(),
     'execution':'Registered wrapper, including registry/hash gates and all four mathematical proof programs',
     'external_neighbor_counts_required':False,'arithmetic_backend':'FLINT/Arb shared; independent proof implementations'}
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('FULL CERTIFICATE PASS',cfg['dimension'],cfg['strict_density_upper_bound'],flush=True)
