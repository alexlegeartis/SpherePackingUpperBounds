from pathlib import Path
from fractions import Fraction as F
import json,itertools,hashlib,numpy as np
from flint import fmpz_mat,fmpq_mat
D=Path(__file__).resolve().parent;P=D.parent/'l80_code_decoder_20261117';L=D.parent/'bachoc_nebe80_reconstruction_20261117';z=np.load(L/'lattice80.npz');G=z['gram'];M=z['twice_root_coordinate_basis'];C=z['binary_code'];g=json.loads((P/'geometry83.json').read_text());Mi=fmpq_mat(M.tolist()).inv();GG=fmpz_mat(G.tolist());sources={k:np.array(json.loads((P/f'case{k}.json').read_text())['r'],dtype=np.int64)for k in [101,108,203,301,302]};case_for={'a':101,'b':101,'c':108,'cross_ac':203,'cross_bc':203};receipts=[]
CW={tuple(np.array(bits)@C%2)for bits in itertools.product(range(2),repeat=5)}
for name,case in case_for.items():
 row=g[name];bp=row['block_permutation'];assert sorted(bp)==list(range(10));assert{tuple(np.array(w)[bp])for w in CW}==CW
 # Rebuild the FULL integral automorphism, not just a single transported vector.
 MP=np.zeros_like(M)
 for j in range(10):MP[:,8*bp[j]:8*bp[j]+8]=M[:,8*j:8*j+8]
 Q=fmpq_mat(MP.tolist())*Mi;assert all(x.denominator==1 for r in Q.tolist()for x in r);U=fmpz_mat([[int(x)for x in r]for r in Q.tolist()]);assert abs(int(U.det()))==1 and U*GG*U.transpose()==GG
 lift=np.array((fmpz_mat([sources[case].tolist()])*U).tolist()[0],dtype=np.int64);assert lift.tolist()==row['lift']and(lift%2).tolist()==row['class'];receipts.append({'name':name,'source_case':case,'determinant':int(U.det()),'block_permutation':bp})
a=np.array(g['a']['class']);b=np.array(g['b']['class']);c=np.array(g['c']['class']);R=np.array([a,b,c]);classes={tuple((np.array(bits)@R)%2):bits for bits in itertools.product(range(2),repeat=3)};assert len(classes)==8
assert np.array_equal((a+c)%2,g['cross_ac']['class'])and np.array_equal((b+c)%2,g['cross_bc']['class']);assert np.array_equal((a+b)%2,sources[301])and np.array_equal((a+b+c)%2,sources[302])
# Fractional parts of all squared distances in each parent half-coset.
mods={k:int(r@G@r)%4 for k,r in sources.items()};assert mods=={101:0,108:0,203:2,301:0,302:0}
Q=np.array([[8,0,3],[0,8,3],[3,3,8]],dtype=np.int64);V=[[F(int(x),4)for x in row]for row in Q];assert int(fmpz_mat(Q.tolist()).det())==368
# V−(1/2)I is positive by Sylvester: principal minors6,36,108.
A=Q-2*np.eye(3,dtype=np.int64);assert[int(fmpz_mat(A[:k,:k].tolist()).det())for k in range(1,4)]==[6,36,108]
minimum={};attain={}
for n in itertools.product(range(-4,5),repeat=3):
 if not any(n):continue
 parity=tuple(x%2 for x in n);q=F(sum(n[i]*int(Q[i,j])*n[j]for i in range(3)for j in range(3)),4)
 if parity not in minimum or q<minimum[parity]:minimum[parity]=q;attain[parity]=n
expected={(0,0,0):F(8),(1,0,0):F(2),(0,1,0):F(2),(0,0,1):F(2),(1,1,0):F(4),(1,0,1):F(5,2),(0,1,1):F(5,2),(1,1,1):F(3)};assert minimum==expected
# Outside this box, ||n||²≥25 and nVn>25/2>8, proving global parity minima.
depths={(1,0,0):F(6),(0,1,0):F(6),(0,0,1):F(6),(1,1,0):F(4),(1,0,1):F(11,2),(0,1,1):F(11,2),(1,1,1):F(5)}
assert all(depths[e]+minimum[e]>=8 for e in depths)
# Independent rational83Gram and block determinant.
Gq=fmpq_mat(G.tolist());Rq=fmpq_mat(R.tolist());Vq=fmpq_mat([[str(x)for x in row]for row in V]);cross=Gq*Rq.transpose()/2;bottom=Rq*Gq*Rq.transpose()/4+Vq
rows=[Gq.tolist()[i]+cross.tolist()[i]for i in range(80)]+[(Rq*Gq/2).tolist()[i]+bottom.tolist()[i]for i in range(3)];Gram=fmpq_mat(rows);assert str(Gram.det())=='23/4'
(D/'gram83.json').write_text(json.dumps([[str(x)for x in row]for row in Gram.tolist()])+'\n')
out={'status':'PASS exact automorphisms, classrelations, vertical parity minima,83Gram determinant; conditional on independently audited parent/depth premises','automorphisms':receipts,'case_norm_mod4':mods,'vertical_gram':[[str(x)for x in row]for row in V],'vertical_minima':{''.join(map(str,k)):str(v)for k,v in minimum.items()},'vertical_attainers':{''.join(map(str,k)):v for k,v in attain.items()},'parent_depths_required':{''.join(map(str,k)):str(v)for k,v in depths.items()},'gram_determinant':'23/4','minimum_squared':8,'center_density_squared':f'{2**85}/23','scope':'No decoder absence replay here. Requires parent minimum8 and cases101,108 depth6;203 depth11/2;301 depth4;302 depth5.'}
(D/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
