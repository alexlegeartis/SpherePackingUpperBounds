"""Replay all exact mathematical gates of this registered witness."""
import json,hashlib,subprocess,sys
from pathlib import Path
from fractions import Fraction
if not __debug__:raise RuntimeError('Assertions must remain enabled')
D=Path(__file__).resolve().parent;R=D.parents[2];cfg=json.loads((D/'config.json').read_text())
for name,want in json.loads((D/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((D/name).read_bytes()).hexdigest()==want,name
package=str(D.relative_to(R));row=next(x for x in json.loads((R/'RESULTS.json').read_text()) if x['dimension']==cfg['dimension'])
if row['package']==package:assert Fraction(row['strict_density_upper_bound'])==Fraction(cfg['strict_density_upper_bound'])
else:
 assert package in row.get('earlier_repository_packages',[])
 assert Fraction(row['strict_density_upper_bound'])<=Fraction(cfg['strict_density_upper_bound'])
executed=[]
for script in ['primary_angular','primary_radial','independent_angular','independent_radial']:
 p=subprocess.run([sys.executable,str(D.parents[1]/'profile_enrichment_20261114'/(script+'.py')),str(D)],capture_output=True,text=True)
 (D/(script+'.log')).write_text(p.stdout+p.stderr);print(p.stdout,flush=True)
 assert p.returncode==0,(script,p.stderr)
 executed.append({'script':script,'returncode':p.returncode})
out={'pass':True,'dimension':cfg['dimension'],'strict_density_upper_bound':cfg['strict_density_upper_bound'],'complete_primary_and_independent_replay':True,'executed':executed,'external_neighbor_counts_required':False,'arithmetic_backend':'Independent implementations sharing FLINT/Arb'}
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('FULL CERTIFICATE PASS',cfg['dimension'],cfg['strict_density_upper_bound'],flush=True)
