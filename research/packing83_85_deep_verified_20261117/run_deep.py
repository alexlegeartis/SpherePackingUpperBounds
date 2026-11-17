from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile,shutil
p=Path(__file__).resolve().parent;ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
if not __debug__:raise RuntimeError('Assertions required')
for f,h in json.load(open(p/'SHA256SUMS.json')).items():assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
with tempfile.TemporaryDirectory()as td:
 w=Path(td)/'work';shutil.copytree(p,w)
 def run(cmd):subprocess.run(list(map(str,cmd)),cwd=w,check=True)
 run([sys.executable,w/'run_all.py']+(['--full']if args.full else[]));run([sys.executable,w/'deep_extension/check.py'])
 for f in ['check_alt.py','check84D4.py','check85D5.py']:
  run([sys.executable,w/'independent_audit'/f])
 if args.full:
  q=w/'deep_extension';run(['g++','-O3',q/'depth4.cpp','-o',q/'depth4'])
  for f in sorted(q.glob('input*.txt')):
   ans=subprocess.run([str(q/'depth4'),str(f),'10000000000'],capture_output=True,text=True,check=True);assert json.loads(ans.stdout)['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT'
 print('PASS83–85 deep-plane extension; '+('full exclusions rerun'if args.full else'quick mode inherits archived exclusions'))
