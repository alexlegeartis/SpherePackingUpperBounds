from pathlib import Path
import subprocess,sys,json,gzip,hashlib,tempfile,shutil,argparse
if not __debug__:raise RuntimeError('Do not disable assertions')
p=argparse.ArgumentParser();p.add_argument('--quick',action='store_true');args=p.parse_args();D=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='lattice36-') as tmp:
 T=Path(tmp);stored=json.loads((D/'STORED_FILES.json').read_text());manifest=json.loads((D/'MANIFEST.json').read_text())
 for f,h in manifest.items():
  if f in stored:
   z=stored[f];raw=(D/z['stored']).read_bytes();assert hashlib.sha256(raw).hexdigest()==z['sha256'];raw=gzip.decompress(raw)
  else:raw=(D/f).read_bytes()
  assert hashlib.sha256(raw).hexdigest()==h;out=T/f;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(raw)
 jobs=[['check.py'],['source37/cap_independent.py']]
 if not args.quick:jobs += [['source37/replay_trace.py','main'],['replay/certify_annular.py'],['check.py'],['replay/finalize.py']]
 for job in jobs:subprocess.run([sys.executable,str(T/job[0]),*job[1:]],check=True,cwd=T)
 print('QUICK EXACT CHECKS PASS; NOT FULL REPLAY' if args.quick else 'COMPLETE36D PORTABLE REPLAY PASS')
