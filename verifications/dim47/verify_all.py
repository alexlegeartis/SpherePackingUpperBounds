"""Replay all mathematical inputs before exact transfer and rigorous rounding."""
import hashlib,json,subprocess,sys
from pathlib import Path
if not __debug__: raise RuntimeError('Do not use Python -O')
D=Path(__file__).parent
expected={
'adaptive_cyclic40_radial111.json':'fac908429ec69fdbbc20b28b1d23fa792822cf0c',
'cyclic_angular40_exact.json':'33e54a15e6fcc9684847efb97235a0cafd91a457',
'independent_cyclic40_angular.py':'f6b909546b5da35790e09df5cd3d97e0b3acb23c',
'independent_cyclic40_radial.py':'37e16d3101cfac350458b6510b5e87d13e8166fb',
'independent_cyclic40_domain.py':'4f7f257ae72a13d90dffe07ce442da1ef5c4afca',
'cyclic40_radial111_exact.json':'8b648bfacdece23e0f1c2a87b0dd8f2299fbc03d'}
for name,want in expected.items():
 data=(D/'continuum_three_point'/name).read_bytes();assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==want
for short in ['angular','radial','domain']:
 print('Replaying source48',short,flush=True)
 p=subprocess.run([sys.executable,'continuum_three_point/independent_cyclic40_'+short+'.py'],cwd=D,capture_output=True,text=True)
 (D/('source_'+short+'_replay.log')).write_text(p.stdout+p.stderr)
 print(p.stdout,flush=True);assert p.returncode==0,(short,p.stderr)
subprocess.run([sys.executable,'verify_transfer.py'],cwd=D,check=True)
(D/'full_replay_result.json').write_text(json.dumps({'pass':True,'source48_full_replay':True,'exact_transfer_replay':True,'dimension':47,'strict_upper':'0.000000616879891','mathematical_inputs':['Cohn-de Laat-Salmon cyclic continuum theorem','Gegenbauer addition theorem','Fourier restriction identity for Schwartz functions']},indent=2))
print('FULL_PROOF_CHAIN_PASS')
