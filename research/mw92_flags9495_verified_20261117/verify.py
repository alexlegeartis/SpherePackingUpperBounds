from pathlib import Path
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required; do not use Python -O')
p=Path(__file__).resolve().parent
for name,h in json.load(open(p/'MANIFEST.json')).items():assert hashlib.sha256((p/name).read_bytes()).hexdigest()==h,name
for script in ['check.py','check95.py','replay_parent.py']:
 subprocess.run([sys.executable,str(p/script)],check=True)
print('PASS 94/95 exact constructions, parent replay and outward rounding')
