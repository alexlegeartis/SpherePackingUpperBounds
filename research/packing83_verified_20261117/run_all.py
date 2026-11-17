from pathlib import Path
import argparse,hashlib,json,shutil,tempfile,subprocess,sys
p=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
if not __debug__:raise RuntimeError('Assertions required')
for fn,h in json.loads((p/'SHA256SUMS.json').read_text()).items():assert hashlib.sha256((p/fn).read_bytes()).hexdigest()==h,fn
with tempfile.TemporaryDirectory(prefix='packing83_')as tmp:
 w=Path(tmp)/'work';shutil.copytree(p,w)
 def run(cmd):
  print('RUN',*map(str,cmd),flush=True);return subprocess.run(list(map(str,cmd)),check=True,cwd=w)
 run([sys.executable,w/'run82.py']+(['--full']if args.full else[]))
 P=w/'l80_code_decoder_20261117';A=w/'l80_case301302_independent_20261117'
 for f in[A/'check_maps.py',w/'l80_83_geometry_independent_20261117/check.py',P/'build83.py',A/'check_density.py']:run([sys.executable,f])
 for ix in[301,302]:
  assert 'EXHAUSTIVE_NO_VECTOR' in(A/f'search{ix}.log').read_text()
  assert json.loads((P/f'case{ix}_run.json').read_text())['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT'
 if args.full:
  for ix,src in[(301,'pair_depth4.cpp'),(302,'pair_depth5.cpp')]:
   run(['g++','-O3','-std=c++17',P/src,'-o',P/f'run{ix}']);q=subprocess.run([str(P/f'run{ix}'),str(P/f'case{ix}.txt'),'10000000000'],capture_output=True,text=True,check=True);print(q.stdout);assert json.loads(q.stdout)['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT'
   run(['g++','-O3','-std=c++17',A/f'search{ix}.cpp','-o',A/f'run{ix}']);q=subprocess.run([str(A/f'run{ix}'),str(A/f'case{ix}.txt')],capture_output=True,text=True,check=True);print(q.stdout);assert 'EXHAUSTIVE_NO_VECTOR' in q.stdout
 print('PASS83 exact target maps, automorphisms, congruences, mixed-depth Gram/density, inherited81/82 proofs and archived completed searches'+('; all primary/independent searches freshly rerun'if args.full else '; quick mode does not rerun exclusions'))
