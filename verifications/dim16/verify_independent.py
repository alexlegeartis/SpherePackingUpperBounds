"""Portable aggregate independent replay, separate from verify_all.py."""
if not __debug__:
    raise RuntimeError('Assertions are required; do not use -O or -OO.')
from pathlib import Path
import os,subprocess,sys,json,shutil
D=Path(__file__).resolve().parent
P=D/'count/dim16_log_count_combined_independent_20261111'
env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
subprocess.run([sys.executable,str(P/'run_all.py')],cwd=D,env=env,check=True)
result=json.loads((P/'COMBINED_RESULT.json').read_text())
assert result['status']=='PASS'
shutil.copyfile(P/'COMBINED_RESULT.json',D/'INDEPENDENT_RESULT.json')
print('PASS: complete independent audit; see INDEPENDENT_RESULT.json',flush=True)
