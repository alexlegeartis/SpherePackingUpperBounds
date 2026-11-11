"""Replay frozen solid-harmonic cyclic certificates, fail closed."""
import subprocess,sys
from pathlib import Path
assert __debug__, 'Do not disable assertions.'
root=Path(__file__).resolve().parent
for n in [17,18,19,20]:
    subprocess.run([sys.executable,str(root/'verify.py'),str(root/f'certificate{n}.json')],check=True)
print('All four primary exact/interval certificates passed.')
