from pathlib import Path
import json,numpy as np
from flint import fmpq_mat as Q,fmpz_mat as Z
p=Path(__file__).resolve().parent;z=np.load(p.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');M=Q(z['twice_root_coordinate_basis'].tolist());Mi=M.inv();G=Q(z['gram'].tolist());K=json.loads((p/'isometries.json').read_text())['isometries'];items=json.loads((p/'deep_plane.json').read_text());classes=[];report=[]
for rec in items:
 U=K[rec['isometry_index']];bp=rec['block_permutation'];assert sorted(bp)==list(range(10));T=Q(80,80)
 for b in range(10):
  for i in range(8):
   for j in range(8):T[8*b+i,8*bp[b]+j]=U[i][j]
 V=M*T*Mi
 assert all(x.denominator==1 for x in V.entries()) and V*G*V.transpose()==G and abs(V.det())==1
 src=json.loads((p.parent/f"l80_code_decoder_20261117/case{rec['source_case']}.json").read_text())['r'];lift=[int(x)for x in (Q([src])*V).entries()]
 assert lift==rec['lift'];cl=[x%2 for x in lift];assert cl==rec['class'] and any(cl);classes.append(cl)
 report.append({'source_case':rec['source_case'],'isometry_index':rec['isometry_index'],'block_permutation':bp,'determinant':str(V.det())})
assert classes[0]!=classes[1] and all((a+b+c)%2==0 for a,b,c in zip(*classes))
(p/'plane_independent_transfer_check.json').write_text(json.dumps({'pass':True,'nonzero_distinct_classes':3,'rank':2,'transfers':report,'premise':'Depth >=6 for source cases101/108 remains inherited, not rerun here.'},indent=2)+'\n');print('PASS THREE EXACT COMPOSITE TRANSFERS; XOR PLANE RANK2')
