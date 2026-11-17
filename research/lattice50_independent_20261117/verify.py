from pathlib import Path
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
for f,h in json.loads((p/'INPUTS_MANIFEST.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
subprocess.run([sys.executable,'-B',str(p/'run_all.py')],check=True)
