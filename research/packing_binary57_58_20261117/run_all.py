from pathlib import Path
import argparse,subprocess,sys,tempfile,shutil,json,hashlib
p=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true',help='Recompute both exhaustive exclusions (about 25 minutes)');args=ap.parse_args()
if not __debug__:raise RuntimeError('Do not run Python with -O')
for f,h in json.loads((p/'SHA256SUMS.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
with tempfile.TemporaryDirectory(prefix='binary57_58_') as tmp:
 t=Path(tmp)
 for d in p.iterdir():
  if d.is_dir():shutil.copytree(d,t/d.name)
 def run(cmd):
  print('RUN',*map(str,cmd),flush=True);subprocess.run(list(map(str,cmd)),check=True,cwd=t)
 primary=t/'cm57_deeper_20261117';ind=t/'cm57_order2_independent_20261117';old=t/'cm57_independent_20261117'
 for f in [old/'input/parent/verify_parent.py',old/'check_setup.py',primary/'verify_binary_setup.py',ind/'check_setup.py']:
  run([sys.executable,f])
 a=json.loads((primary/'order2_527_complete.json').read_text());b=json.loads((ind/'search.log').read_text().strip().splitlines()[-1])
 assert a['status']==b['status']=='EXHAUSTIVE_NO_VECTOR' and a['nodes']<10000000000 and b['nodes']<6000000000
 flags=['g++','-O3','-std=c++17','-ffp-contract=off','-fno-fast-math','-frounding-math']
 run(flags+[ind/'controls_search.cpp','-o',ind/'controls_search']);run([sys.executable,ind/'test_controls.py'])
 if args.full:
  run(flags+[primary/'verify_order2.cpp','-o',primary/'verify_order2'])
  run(flags+[ind/'independent_search.cpp','-o',ind/'independent_search'])
  for name,cmd in [('primary',[primary/'verify_order2',primary/'order2_527_bits.txt','10000000000','fast','bits']),('independent',[ind/'independent_search',ind/'input.txt'])]:
   print('EXHAUSTIVE',name,flush=True)
   q=subprocess.run(list(map(str,cmd)),capture_output=True,text=True,check=True);print(q.stdout,flush=True)
   result=json.loads(q.stdout.strip().splitlines()[-1]);assert result['status']=='EXHAUSTIVE_NO_VECTOR'
 print('PASS: exact setup, parent provenance, recorded completed proofs and 40 independent controls'+('; both exclusions freshly recomputed' if args.full else '; exclusions not freshly recomputed in quick mode'))
