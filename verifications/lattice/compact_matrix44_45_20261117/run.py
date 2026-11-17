import argparse,tempfile,shutil,subprocess,sys,hashlib,json
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['quick','full'],default='quick');args=p.parse_args();D=Path(__file__).resolve().parent
for name,h in json.load(open(D/'MANIFEST.json')).items():assert hashlib.sha256((D/name).read_bytes()).hexdigest()==h,name
with tempfile.TemporaryDirectory()as td:
 T=Path(td)/'proof';shutil.copytree(D,T)
 for f in T.rglob('*.py'):
  s=f.read_text().replace('research/lattice4445_independent_20261117/input44',str(T/'input44')).replace('research/lattice4445_independent_20261117/input45',str(T/'input45'));f.write_text(s)
 def run(name,*args):subprocess.run([sys.executable,str(T/name),*args],check=True,cwd=T)
 run('input45/restrict44_45.py');run('audit_links.py');run('source46/check_exact.py')
 for n in [44,45]:
  run(f'input{n}/check_fourier'+('45'if n==45 else'')+'.py')
 run('input44/certify_cap_exact.py')
 if args.mode=='full':
  run('replay44.py','main')
  for n in [44,45]:
   run(f'input{n}/root_repair{n}/check_tail.py');run(f'input{n}/root_repair{n}/certify_annular.py')
  run('input44/normalize.py');run('input45/normalize45.py')
  run('source46/input/root_repair46/check_tail.py');run('source46/replay_trace.py','main');run('source46/scalar_replay.py');run('source46/check_exact.py')
 print('PASS',args.mode)
