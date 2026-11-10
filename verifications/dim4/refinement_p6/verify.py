"""Replay the frozen 4D certificate in an isolated temporary directory."""
from pathlib import Path
from fractions import Fraction
import json,subprocess,sys,tempfile,shutil
assert __debug__, 'Do not disable assertions'
root=Path(__file__).resolve().parent
scripts=['certify_radial.py','certify_psd.py','certify_count_pair_arb.py','certify_far_bernstein.py','certify_cyclic_bernstein.py','certify_density.py']
with tempfile.TemporaryDirectory(prefix='sphere4-') as d:
 d=Path(d)
 shutil.copyfile(root/'witness.json',d/'witness.json')
 for name in scripts:shutil.copyfile(root/name,d/name)
 for name in scripts:
  arg='witness_psd_exact.json' if name in scripts[2:5] else 'witness.json'
  subprocess.run([sys.executable,name,arg],cwd=d,check=True)
 def read(suffix):return json.loads((d/('witness'+suffix+'.json')).read_text())
 a=read('');audited=json.loads((root/'audited_snapshot.json').read_text());assert all(a[k]==audited[k] for k in a)
 assert (a['R'],a['J'],a['p'],a['K'])==(1.665,1,6,16)
 assert a['count_radii']==[1.008096,1.017445,1.025977,1.03695]
 assert a['count_bounds']==[24.0,25.0,26.0,27.0]
 assert all(Fraction(str(w))>=0 for w in a['count_weights'])
 psd=read('_psd_exact');assert all(x['positive_definite'] for x in psd['checks'])
 assert psd['matrices']==json.loads((root/'witness_psd_exact.json').read_text())['matrices']
 rad=read('_radial_exact');assert rad['fourier_positive_roots']==0 and rad['fourier_origin_positive']
 assert rad['tail_roots_after_conservative_threshold']==0 and rad['tail_start_negative']
 for suffix in ['_count_pair_arb','_far_bernstein','_cyclic_bernstein']:assert read(suffix)['verified'],suffix
 assert read('_density_exact')['below_0.633430327']
print('PASS: dimension 4 density < 0.633430327, using the cited external neighbor-count theorem')
