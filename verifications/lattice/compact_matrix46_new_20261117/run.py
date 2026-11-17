import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
if not __debug__:raise RuntimeError('Do not use -O')
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['quick','full'],default='quick');args=p.parse_args();D=Path(__file__).resolve().parent
for rel,h in json.load(open(D/'MANIFEST.json')).items():assert hashlib.sha256((D/rel).read_bytes()).hexdigest()==h,rel
with tempfile.TemporaryDirectory(prefix='lattice46-proof-')as t:
 W=Path(t)/'proof';shutil.copytree(D/'proof',W)
 def run(*cmd):subprocess.run([sys.executable,*map(str,cmd)],check=True,cwd=W)
 run(W/'input/certify_cap_exact.py');run(W/'check_exact.py')
 if args.mode=='full':
  run(W/'input/root_repair46/check_tail.py');run(W/'scalar_replay.py');run(W/'check_exact.py');run(W/'replay_trace.py','main')
 print('PASS',args.mode,'46D density <= 7.766101459e-7')
 if args.mode=='quick':print('Quick mode checks hashes and exact algebra only. Use --mode full to replay all numerical interval enclosures.')
