from pathlib import Path
import subprocess,sys,json,hashlib
D=Path(__file__).resolve().parent
for path,digest in json.loads((D/'INPUT_SHA256.json').read_text()).items():
 assert hashlib.sha256((D/path).read_bytes()).hexdigest()==digest,path
for command in [['check_exact.py'],['geometry/verify_caps.py','51'],['geometry/verify_edges.py','51'],['scalar/certify_annuli.py'],['check_exact.py'],['replay_geometry.py','51']]:
 subprocess.run([sys.executable,*command],cwd=D,check=True)
a=json.loads((D/'main_independent_result.json').read_text()); assert a['nodes']==106453 and a['leaves']==53227
x=json.loads((D/'exact_result.json').read_text());assert x['direct_tail_shift_all_nonpositive']
assert a['input_sha']==x['frozen_sha256']
print('PASS unrestricted 51D packing density < 1.37012082e-7; complete exact/interval replay')
