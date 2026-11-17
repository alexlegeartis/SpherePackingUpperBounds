from pathlib import Path
import sys,tempfile,shutil,json,hashlib,subprocess
if not __debug__:raise RuntimeError('Assertions required')
P=Path(__file__).resolve().parent;G=P.parent/'compact_matrix48_20261117'
for f,h in json.loads((P/'SHA256SUMS.json').read_text()).items():assert hashlib.sha256((P/f).read_bytes()).hexdigest()==h,f
assert G.is_dir(),'The sibling compact_matrix48_20261117 geometry package is required'
with tempfile.TemporaryDirectory()as temp:
 D=Path(temp)/'proof';shutil.copytree(P,D);subprocess.run([sys.executable,'-B',str(D/'full_replay.py'),'--geometry',str(G),*sys.argv[1:]],check=True)
