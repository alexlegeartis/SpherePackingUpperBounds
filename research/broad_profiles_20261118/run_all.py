from pathlib import Path
import subprocess,sys
from concurrent.futures import ThreadPoolExecutor
R=Path(__file__).resolve().parents[2]
def run(n):subprocess.run([sys.executable,str(R/f'verifications/dim{n}/refinement_broad20261118/verify.py')],check=True)
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,range(17,24)))
print('PASS: all seven registered broader-profile certificates')
