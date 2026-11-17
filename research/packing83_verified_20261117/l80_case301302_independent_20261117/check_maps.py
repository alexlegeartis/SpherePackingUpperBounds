from pathlib import Path
import json,numpy as np,hashlib
from flint import fmpz_mat as Z,fmpq_mat as Q
D=Path(__file__).resolve().parent;P=D.parent/'l80_code_decoder_20261117';L=D.parent/'bachoc_nebe80_reconstruction_20261117';d=np.load(L/'lattice80.npz');G=Z(d['gram'].tolist());M=Z(d['twice_root_coordinate_basis'].tolist());E=d['e8_twice_euclidean_basis']
old=json.loads((D.parent/'l80_case101_independent_20261117/case101.json').read_text());oldtxt=list(map(int,(D.parent/'l80_case101_independent_20261117/input.txt').read_text().split()))
# Layout pinned and fully reconstructed in the previous independent mapping proof.
start=1+9121*8+9121+64+256
out={}
for n,mod,cut in [(301,0,384),(302,0,512)]:
 c=json.loads((D/f'case{n}.json').read_text());r=Z([c['r']]);assert len(c['r'])==80 and set(c['r'])<={0,1};assert G==Z(c['gram']);assert int((r*G*r.transpose())[0,0])%4==mod
 vals=list(map(int,(D/f'case{n}.txt').read_text().split()));assert len(vals)==len(oldtxt);assert vals[:start]==oldtxt[:start] and vals[start+80:]==oldtxt[start+80:]
 X=np.array((r*M).tolist(),dtype=object).reshape(10,8)@E;assert vals[start:start+80]==list(map(int,X.flat));assert np.array_equal(X,c['target_euclidean_times8'])
 out[str(n)]={'target_norm':str((r*G*r.transpose())[0,0]),'norm_mod4':mod,'cutoff_numerator':cut,'input_sha256':hashlib.sha256((D/f'case{n}.txt').read_bytes()).hexdigest()}
out['status']='PASS new exact targets, common decoder inputs and integer-distance congruences';(D/'mapping_result.json').write_text(json.dumps(out,indent=2));print(out)
