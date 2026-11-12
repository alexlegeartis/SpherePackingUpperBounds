#!/usr/bin/env python3
"""Parse archived catalogue Grams and prepare exact congruence reductions."""
import re,json,sys,pathlib,hashlib
import numpy as np
from fpylll import IntegerMatrix,GSO,LLL
HERE=pathlib.Path(__file__).resolve().parent

def sec(s,k):
 m=re.search(r'<a\s+NAME="'+k+r'[^\"]*">.*?</a>(.*?)(?=<p><li>|<a\s+NAME=|<hr>|$)',s,re.S|re.I)
 return re.sub(r'<[^>]+>',' ',m.group(1)).strip() if m else ''

def parse(name):
 s=(HERE/'sources'/f'{name}.html').read_text()
 n=int(sec(s,'DIMENSION'))
 values=list(map(int,re.findall(r'-?\d+',sec(s,'GRAM'))))
 if len(values)==n*n+2:
  assert values[:2]==[n,n];G=np.array(values[2:],dtype=np.int64).reshape(n,n)
 elif len(values)==n*n+1:
  assert values[0]==n;G=np.array(values[1:],dtype=np.int64).reshape(n,n)
 elif len(values)==n*(n+1)//2+2:
  assert values[:2]==[n,0];G=np.zeros((n,n),np.int64);ii,jj=np.tril_indices(n);G[ii,jj]=values[2:];G[jj,ii]=values[2:]
 else:raise ValueError((name,n,len(values),values[:4]))
 assert np.array_equal(G,G.T)
 assert min(np.linalg.eigvalsh(G))>0
 return G,{k:re.sub(r'\s+',' ',sec(s,k)) for k in ['DIMENSION','DET','DIVISORS','MINIMAL_NORM','KISSING_NUMBER']}

def reduce_gram(G):
 n=len(G);A=IntegerMatrix.from_matrix(G.tolist());U=IntegerMatrix.identity(n)
 M=GSO.Mat(A,U=U,gram=True,flags=GSO.INT_GRAM,float_type='long double',update=True)
 LLL.Reduction(M,delta=.99)()
 H=np.array([[int(A[i,j]) for j in range(n)]for i in range(n)],np.int64)
 H=np.tril(H)+np.tril(H,-1).T
 U=np.array([[int(U[i,j]) for j in range(n)]for i in range(n)],np.int64)
 assert np.array_equal(U@G@U.T,H)
 return H,U

if __name__=='__main__':
 for name in sys.argv[1:]:
  G,meta=parse(name);H,U=reduce_gram(G);out=HERE/name;out.mkdir(exist_ok=True)
  np.savez(out/'seed.npz',G=G,H=H,U=U)
  np.savetxt(out/'Q.txt',H,fmt='%d');np.savetxt(out/'R.txt',np.linalg.cholesky(H).T,fmt='%.17g')
  meta.update({'source_sha256':hashlib.sha256((HERE/'sources'/f'{name}.html').read_bytes()).hexdigest(),'logdet':float(np.linalg.slogdet(G)[1]),'reduced_diagonal':H.diagonal().tolist(),'GSO_diagonal':(np.linalg.cholesky(H).diagonal()**2).tolist()})
  (out/'metadata.json').write_text(json.dumps(meta,indent=2));print(name,json.dumps(meta),flush=True)
