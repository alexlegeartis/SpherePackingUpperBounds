from pathlib import Path
import json,itertools,hashlib,numpy as np
from fractions import Fraction as F
from flint import fmpq_mat as Q
D=Path(__file__).resolve().parent;P=D.parent/'deep_extension';L=D.parent/'bachoc_nebe80_reconstruction_20261117';z=np.load(L/'lattice80.npz');G=Q(z['gram'].tolist());M=Q(z['twice_root_coordinate_basis'].tolist());E=z['e8_twice_euclidean_basis'];Mi=M.inv();c=json.loads((P/'candidate84.json').read_text());iso=json.loads((D.parent/'l80_plane82_independent_20261117/isometries.json').read_text())['isometries'];cl=[]
for it in c['generators']:
 src=json.loads((D.parent/'l80_code_decoder_20261117'/f"case{it['source_case']}.json").read_text());assert Q(src['gram'])==G
 U=Q(iso[it['isometry_index']]);rows=[]
 for row in M.tolist():
  blocks=[None]*10
  for j,b in enumerate(it['block_permutation']):blocks[b]=(Q([row[8*j:8*j+8]])*U).tolist()[0]
  rows.append(sum(blocks,[]))
 A=Q(rows)*Mi;assert all(x.denominator==1 for x in A.entries());assert abs(A.det())==1 and A*G*A.transpose()==G
 v=Q([src['r']])*A;assert v==Q([it['lift']]);assert [int(x)%2 for x in v.entries()]==it['class'];cl.append(it['class'])
assert hashlib.sha256((D.parent/'l80_code_decoder_20261117/case101.json').read_bytes()).hexdigest()=='bced37fa26cf53923af5f5fc59379203939f3b96f3692e8945110c6c2351022b'
assert hashlib.sha256((D.parent/'l80_code_decoder_20261117/case101.txt').read_bytes()).hexdigest()=='5cb83cfe951a5b4968aa6e32d817fa2952c3a1668017fb5d7a3e0f4bcda76d32'
assert hashlib.sha256((L/'lattice80.npz').read_bytes()).hexdigest()=='e157bdeba01155056be8b9f844550353b5bf9a34340649eb62d1ef62ae69ba81'
# Check advertised depth transfers for all five known classes against exact stored orbit records.
records=list(c['known'].values())
# Each needed match independently verifies full80 isometry, not family recognition.
for mask in [1,2,3,4,5,6,8,9,10,12]:
 target=[sum(cl[j][i]for j in range(4)if mask>>j&1)%2 for i in range(80)]
 rec=next(r for r in records if r['class']==target);src=json.loads((D.parent/'l80_code_decoder_20261117'/f"case{rec['source_case']}.json").read_text());U=Q(iso[rec['isometry_index']]);rows=[]
 for row in M.tolist():
  blocks=[None]*10
  for j,b in enumerate(rec['block_permutation']):blocks[b]=(Q([row[8*j:8*j+8]])*U).tolist()[0]
  rows.append(sum(blocks,[]))
 A=Q(rows)*Mi;assert all(x.denominator==1 for x in A.entries())and abs(A.det())==1 and A*G*A.transpose()==G
 assert [int(x)%2 for x in (Q([src['r']])*A).entries()]==target
 print('mask',mask,'source',rec['source_case'])
V=Q([[2 if i==j else 1 for j in range(4)]for i in range(4)]);assert V.det()==5
W=V-Q([["1/2" if i==j else 0 for j in range(4)]for i in range(4)]);assert all(Q([[W[i,j]for j in range(k)]for i in range(k)]).det()>0 for k in range(1,5))
mins=[None]*16
for a in itertools.product(range(-4,5),repeat=4):
 if not any(a):continue
 val=sum(F(int(V[i,j].numerator),int(V[i,j].denominator))*a[i]*a[j]for i in range(4)for j in range(4));mask=sum((a[i]%2)<<i for i in range(4));mins[mask]=val if mins[mask]is None else min(mins[mask],val)
depths=[0]+[6 if m.bit_count()<=2 else 4 for m in range(1,16)];assert all(a+b>=8 for a,b in zip(mins,depths));assert 25*F(1,2)>8
# all new parity targets bind to exact Gram, residue inputs equal parent decoder except target.
base=list(map(int,(D.parent/'l80_code_decoder_20261117/case101.txt').read_text().split()));n=base[0];start=1+9*n+64+256
for mask in [7,11,13,14,15]:
 r=np.array([sum(cl[j][i]for j in range(4)if mask>>j&1)%2 for i in range(80)]);inp=list(map(int,(P/f'input84_{mask}.txt').read_text().split()));assert inp[:start]==base[:start]and inp[start+80:]==base[start+80:]
 assert inp[start:start+80]==((r@z['twice_root_coordinate_basis']).reshape(10,8)@E).ravel().tolist();assert int(r@z['gram']@r)%4==0
# Full enumeration of local E8 shell by recursive square budget, distinct from product generator.
shell=set()
def rec(prefix,budget,par):
 if len(prefix)==8:
  if sum(prefix)%4==0:shell.add(tuple(prefix))
  return
 for v in range(-4,5):
  if v%2==par and v*v<=budget:rec(prefix+[v],budget-v*v,par)
for par in [0,1]:rec([],24,par)
DS={tuple(base[1+8*i:1+8*i+8])for i in range(n)};assert shell==DS and len(shell)==9121
out={'status':'PASS exact transfer, new target binding, local shell completeness, vertical geometry; global decoders separately pending','vertical_minima':list(map(str,mins)),'determinant':'5','candidate_sha256':hashlib.sha256((P/'candidate84.json').read_bytes()).hexdigest()};(D/'check84_result.json').write_text(json.dumps(out,indent=2));print(out)
