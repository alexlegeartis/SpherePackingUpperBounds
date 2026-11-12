"""Replay every registered proof package. Never run with python -O."""
import subprocess,sys,argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
if not __debug__: raise RuntimeError('Assertions are required')
r=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--jobs',type=int,default=1,help='Parallel dimension chains; each dimension stays sequential (1–8)');args=parser.parse_args()
if not 1<=args.jobs<=8:parser.error('--jobs must be between1 and8')
subprocess.run([sys.executable,str(r/'comparison.py')],check=True)
jobs=[('dim16_reconstruction_minB','verify_all.py')]+[('dim10','verify_all.py')]+[('dim14','verify.py')]+[('dim13','verify.py')]+[('dim12','verify.py')]+[('dim11','verify.py')]+[('dim15','verify.py')]+[('dim46','refinement_extended_pair/verify.py')]+[('dim47','refinement_degree47/verify.py')]+[('dim46','refinement_degree47/verify.py')]+[('dim47','refinement_joint/verify.py')]+[('dim46','refinement_joint/verify.py')]+[('dim47','refinement_angular/verify.py')]+[('dim46','refinement_angular/verify.py')]+[('dim21','refinement_gaussian/verify.py')]+[('dim29', 'verify.py'), ('dim30', 'verify.py'), ('dim31', 'verify.py'), ('dim33', 'verify.py'), ('dim35', 'verify.py'), ('dim37', 'verify.py')]+[('dim22','verify.py'),('dim23','verify.py'),('dim26','verify.py'),('dim27','verify.py')]+[('dim25', 'verify.py'), ('dim28', 'verify.py'), ('dim32', 'verify.py'), ('dim34', 'verify.py'), ('dim36', 'verify.py'), ('dim38', 'verify.py'), ('dim39', 'verify.py'), ('dim40', 'verify.py')]+[('dim41','verify.py'),('dim42','verify.py'),('dim43','verify.py'),('dim44','verify.py'),('dim6','refinement_patch5/verify.py'),('dim6','refinement_mixture99/verify.py'),('dim6','verify.py'),('dim45','verify.py'),('dim46','verify.py'),('dim47','verify.py'),('dim17','refinement_repaired/verify.py'),('dim19','refinement_cyclic/verify.py'),('dim20','refinement_cyclic/verify.py'),('dim17','verify.py'),('dim18','refinement_cyclic/verify.py'),('dim18','verify.py'),('dim5','refinement_p4/verify.py'),('dim4','refinement_p6/verify.py'),('dim5','refinement_pd2/verify.py'),('dim4','refinement_p4/verify.py'),('dim5','verify.py'),('dim4','verify.py'),('dim21','verify.py'),('dim21','independent_verify.py'),('dim19','verify.py'),('dim19','independent_verify.py'),('dim20','verify.py'),('dim20','independent_verify.py')]+[('dim48','continuum_three_point/independent_cyclic40_'+x+'.py') for x in ['angular','domain','radial']]
jobs.append(('dim7','verify.py'))
jobs.append(('dim5','refinement_guarded/verify.py'))
jobs.append(('dim5','refinement_lens/verify.py'))
jobs.append(('dim5','refinement_spherical_lens/verify.py'))
def replay(chain):
 for directory,script in chain:
  print(directory,script,flush=True)
  subprocess.run([sys.executable,script],cwd=r/'verifications'/directory,check=True)
if args.jobs==1:
 replay(jobs)
else:
 groups={}
 for directory,script in jobs:groups.setdefault(directory,[]).append((directory,script))
 with ThreadPoolExecutor(max_workers=args.jobs) as pool:
  futures=[pool.submit(replay,chain) for chain in groups.values()]
  for future in as_completed(futures):future.result()
from fractions import Fraction
import json
a=json.loads((r/'verifications/dim48/continuum_three_point/independent_cyclic40_radial.json').read_text())
assert Fraction(a['density_upper']) < Fraction('4.13385739314e-7')
subprocess.run([sys.executable,str(r/'verifications/origin_slack/verify_algebra.py')],check=True)
subprocess.run([sys.executable,str(r/'verifications/cyclic_origin_slack/verify_algebra.py')],check=True)
print('PASS: all registered proof replays and strict decimal rounding')
