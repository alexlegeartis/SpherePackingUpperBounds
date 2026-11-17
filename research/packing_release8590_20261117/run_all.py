from pathlib import Path
import subprocess,sys,json
if not __debug__:raise RuntimeError('Run without -O.')
p=Path(__file__).resolve().parent
for row in json.load(open(p/'RESULTS.json'))['rows']:
 subprocess.run([sys.executable,str(p/'packing_resume_20261117/verify_candidate.py'),row['file']],check=True)
for f in ['verify182.py','verify192_85.py','verify200_89.py']:
 subprocess.run([sys.executable,str(p/'audit208wide_20261117'/f)],check=True)
print('PASS all six primary and independent checks')
