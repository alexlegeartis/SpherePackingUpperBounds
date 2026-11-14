"""Exact 34D target: five cosets of a size-4^7 distance-eight subcode.

Parent codes are extended beyond an odd-like cyclic code by coordinate moves.
No floating arithmetic is used in code or clique verification.
"""
from pathlib import Path
import json,itertools,time
import numpy as np
P=Path(__file__).resolve().parent
mul4=np.array([[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]],dtype=np.uint8)
# Polynomial generator, low degree first; its eight roots are the 4-cyclotomic
# orbits {1,4,13,16} and {3,5,12,14} among the nonzero residues modulo 17.
g=np.array([1,2,0,2,2,2,0,2,1],dtype=np.uint8)
even=np.r_[g,0]^np.r_[0,g]
G9=np.array([np.pad(g,(i,8-i))for i in range(9)],dtype=np.uint8)
G8=np.array([np.pad(even,(i,7-i))for i in range(8)],dtype=np.uint8)
def words(G):
 C=np.zeros((1,17),dtype=np.uint8)
 for v in G:C=np.concatenate([C^mul4[z,v]for z in range(4)])
 return C
assert np.count_nonzero(words(G8)[1:],axis=1).min()==8
assert np.count_nonzero(words(G9)[1:],axis=1).min()==7
def find_clique(adj,target):
 def go(chosen,available):
  if len(chosen)>=target:return chosen
  if len(chosen)+available.bit_count()<target:return None
  while available:
   bit=available&-available;v=bit.bit_length()-1;available^=bit
   answer=go(chosen+[v],available&adj[v])
   if answer:return answer
  return None
 return go([0],adj[0])
start=time.time();records=[];rng=np.random.default_rng(20261114)
moves=[]
for coordinate in range(17):
 e=np.zeros(17,dtype=np.uint8);e[coordinate]=1;moves.append(e)
for i,j in itertools.combinations(range(17),2):
 for value in (1,2,3):
  e=np.zeros(17,dtype=np.uint8);e[i]=1;e[j]=value;moves.append(e)
for _ in range(256):
 ids=np.sort(rng.choice(17,3,replace=False));e=np.zeros(17,dtype=np.uint8);e[ids]=[1,int(rng.integers(1,4)),int(rng.integers(1,4))];moves.append(e)
for trial in range(4):
 # Random invertible row operations give different codimension-one subcodes.
 T=G8.copy()
 if trial:
  for _ in range(30):
   a,b=rng.choice(8,2,replace=False);T[a]^=mul4[int(rng.integers(1,4)),T[b]]
 S=T[:7];C=words(S)
 for coordinate,e in enumerate(moves):
  extensions=[T[7],G9[0],e];mins=[]
  for index in range(64):
   a=index;v=np.zeros(17,dtype=np.uint8)
   for r in extensions:v^=mul4[a%4,r];a//=4
   w=np.count_nonzero(C^v,axis=1);mins.append(int(w.min()))
  good=[i for i in range(1,64)if mins[i]>=8]
  adj=[sum(1<<(v^a)for a in good)for v in range(64)]
  clique=find_clique(adj,5)
  records.append({'trial':trial,'move_index':coordinate,'move':e.tolist(),'compatible_differences':len(good),'five_cosets':clique})
  if clique:
   np.savez(P/'new34_code_candidate.npz',S=S,extensions=np.array(extensions),clique=np.array(clique),minima=np.array(mins))
   print('FOUND',trial,coordinate,clique,flush=True)
   (P/'cyclic17_cosets_results.json').write_text(json.dumps(records,indent=2)+'\n');raise SystemExit(0)
 print('trial',trial,'parents',len(moves),'reject five cosets','seconds',time.time()-start,flush=True)
(P/'cyclic17_cosets_results.json').write_text(json.dumps(records,indent=2)+'\n')
