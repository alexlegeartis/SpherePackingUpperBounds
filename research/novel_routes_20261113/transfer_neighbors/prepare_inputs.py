"""Recreate catalogue inputs and exact BKZ transformation witnesses.

Fetch is optional: the primary source is archived in sources/2_dim32.dat.
Run with fpylll/cysignals available. Reduction is discovery, while the saved
integral congruences and unimodularity are independently checked exactly.
"""
from pathlib import Path
import re,hashlib,json
import numpy as np,sympy as sp
from fpylll import IntegerMatrix,GSO,LLL,BKZ
P=Path(__file__).resolve().parent
source=P/'sources/2_dim32.dat';s=source.read_text();rows=[]
for i,b in enumerate(s.split('Inner Product Matrix:')[1:]):
 a=[]
 for line in b.splitlines():
  if line.startswith('['):
   row=list(map(int,re.findall(r'-?\d+',line)))
   if len(row)==32:a.append(row)
   if len(a)==32:break
 G=np.array(a,dtype=np.int64);assert G.shape==(32,32);assert np.array_equal(G,G.T)
 assert sp.Matrix(G).det()==65536
 Q=IntegerMatrix.from_matrix(G.tolist());U=IntegerMatrix.identity(32)
 M=GSO.Mat(Q,gram=True,U=U,update=True);lll=LLL.Reduction(M);lll();BKZ.Reduction(M,lll,BKZ.Param(block_size=22,max_loops=2))()
 H=np.array([[Q[j,k]for k in range(32)]for j in range(32)],dtype=np.int64)
 T=np.array([[U[j,k]for k in range(32)]for j in range(32)],dtype=np.int64)
 assert np.array_equal(T@G@T.T,H);assert abs(sp.Matrix(T).det())==1
 for suffix,gram in [('',G),('_bkz',H)]:
  np.savetxt(P/f'Q32_{i}{suffix}.txt',gram,fmt='%d');np.savetxt(P/f'Q32_{i}{suffix}_R.txt',np.linalg.cholesky(gram).T,fmt='%.18g')
 np.savez(P/f'Q32_{i}_bkz.npz',G=G,H=H,T=T)
 rows.append({'index':i,'modular_flag_from_catalogue':re.search(r'Modular\s*=\s*(true|false)',b).group(1)=='true','dimension':32,'determinant':65536,'minimum_source':6,'basis_congruence_exact':True,'basis_determinant':int(sp.Matrix(T).det())})
(P/'source_metadata.json').write_text(json.dumps({'url':'https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/ExtLatDat/2_dim32.dat','sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'lattices':rows,'scope':'The minimum6 is inherited from the primary catalogue, not re-proved by this script.'},indent=2))
