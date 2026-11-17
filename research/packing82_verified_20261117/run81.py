from pathlib import Path
import sys,subprocess,tempfile,shutil,json,hashlib,argparse
p=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
if not __debug__:raise RuntimeError('Assertions must be enabled')
for fn,h in json.loads((p/'SHA256SUMS.json').read_text()).items():assert hashlib.sha256((p/fn).read_bytes()).hexdigest()==h,fn
with tempfile.TemporaryDirectory(prefix='packing81_')as tmp:
 t=Path(tmp)
 for d in p.iterdir():
  if d.is_dir():shutil.copytree(d,t/d.name)
 def run(cmd):
  print('RUN',*map(str,cmd),flush=True);return subprocess.run(list(map(str,cmd)),check=True,cwd=t)
 parent=t/'bachoc_nebe80_reconstruction_20261117';primary=t/'l80_code_decoder_20261117';mapping=t/'l80_case101_mapping_audit_20261117';audit=t/'l80_case101_independent_20261117'
 for f in [parent/'verify_frozen.py',mapping/'check.py',audit/'check_local.py',primary/'build81.py',audit/'check_extension.py']:run([sys.executable,f])
 assert json.loads((primary/'case101_replay.json').read_text())['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT'
 assert 'EXHAUSTIVE_NO_VECTOR independent_reverse_rows nodes=5877535660 leaves=100' in(audit/'search.log').read_text()
 assert json.loads((audit/'result.json').read_text())['status']=='COMPLETE INDEPENDENT AUDIT PASS'
 if args.full:
  for src,exe in [(primary/'decode.cpp',primary/'decode'),(audit/'search.cpp',audit/'search')]:run(['g++','-O3','-std=c++17',src,'-o',exe])
  q=subprocess.run([str(primary/'decode'),str(primary/'case101.txt'),'10000000000'],capture_output=True,text=True,check=True);print(q.stdout);assert json.loads(q.stdout)['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT'
  q=subprocess.run([str(audit/'search'),str(audit/'input.txt')],capture_output=True,text=True,check=True);print(q.stdout);assert 'EXHAUSTIVE_NO_VECTOR independent_reverse_rows nodes=5877535660 leaves=100' in q.stdout
  for j in range(3):
   q=subprocess.run([str(audit/'search'),str(audit/f'control{j}.txt')],capture_output=True,text=True);assert q.returncode==2 and 'COUNTEREXAMPLE' in q.stdout
 print('PASS exact reconstruction, maps, local coverage, Gram/density and archived completed logs'+('; both complete exclusions and3positive controls freshly replayed'if args.full else '; quick mode does not rerun exclusions'))
