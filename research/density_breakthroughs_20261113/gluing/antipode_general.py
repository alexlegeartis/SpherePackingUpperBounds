"""Joint normal-subspace / finite dual-cluster discovery inP48n.
Exact witnesses, saturation and exact density ratios checked on record candidates.
Search and short-dual enumeration are not completeness claims.
"""
from pathlib import Path
from collections import OrderedDict
from fractions import Fraction
import numpy as np,sympy as sp,json,math,time,argparse,re
from sympy.matrices.normalforms import smith_normal_form
from fpylll import Enumeration,GSO,IntegerMatrix,EnumerationError
P=Path(__file__).resolve().parent;ROOT=P.parents[2];ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=500);ap.add_argument('--seed',type=int,default=390038);ap.add_argument('--pool',default='P48n_extended_normal_pool.npz');ap.add_argument('--ranks',default='5,6,7');ap.add_argument('--tag',default='antipode_lower_rank');a=ap.parse_args()
j=np.load(P/a.pool);G=j['G'];baseV=j['V'].astype(np.int16);wanted=[int(k) for k in a.ranks.split(',')];seedrows={6:[j['E8'][:6]],7:[j['E8'][:7]],8:[j['E8']]}
def matrix_after(text,key):
 p=text.index(key);l=text.index('[[',p);r=text.index(']]',l)+2
 return np.array([[int(z)for z in re.findall(r'-?\d+',row)]for row in text[l:r].splitlines()],dtype=np.int16)
EXT=ROOT/'research/density_breakthroughs_20261113/seed_deformations/external_sphere_packing/packings_out/records/dim43-delta226'
seedrows[5]=[]
for file in sorted(EXT.glob('P48n-*.cert.txt')):
 text=file.read_text();GG=matrix_after(text,'ambient Gram G');assert np.array_equal(GG,G)
 seedrows[5].append(matrix_after(text,'coefficient rows w.r.t.'))
for k in [9,10]:seedrows[k]=[np.load(P/f'P48n_best_rank{k}.npz')['B']]
extra=np.vstack([B for k in wanted for B in seedrows[k]]).astype(np.int16);V=np.vstack([baseV,extra]);V*=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)])[:,None];V=np.unique(V,axis=0);N=len(V);norms=np.einsum('bi,ij,bj->b',V.astype(np.int64),G,V.astype(np.int64));lookup={v.tobytes():i for i,v in enumerate(V)};VG=(V.astype(np.int64)@G).astype(float);cache=OrderedDict();rng=np.random.default_rng(a.seed)
target2={5:Fraction(137**43,2**223*3**42),6:Fraction(3**35,2**42),7:Fraction(3**34,2**42),8:Fraction(3**32,2**40),9:Fraction(3**32,2**40*5),10:Fraction(3**31,2**38*25)}
targets={k:float(Fraction(3,2)**(48-k)/target2[k])for k in wanted};start=time.time();record={k:1. for k in wanted};history=[];evaluations=0

def col(i):
 if i in cache:z=cache.pop(i);cache[i]=z;return z
 z=np.rint(VG@V[i].astype(float)).astype(np.int8);assert int(z[i])==int(norms[i]);cache[i]=z
 if len(cache)>256:cache.popitem(last=False)
 return z

def scores(ids):
 if not ids:return norms.astype(float)
 C=np.column_stack([col(i)for i in ids]).astype(float);K=C[ids];q=norms-np.einsum('ij,ij->i',C@np.linalg.inv(K),C);q[q<1e-7]=np.inf;return q

def evaluate(ids):
 global evaluations
 k=len(ids);dim=48-k;B=V[ids].astype(np.int64);K=B@G@B.T;det=int(round(np.linalg.det(K)));J=np.rint(np.linalg.inv(K)*det).astype(np.int64)
 if not np.array_equal(J@K,det*np.eye(k,dtype=np.int64)):return 0.,None
 M=GSO.Mat(IntegerMatrix.from_matrix(J.tolist()),gram=True,update=True);E=Enumeration(M,nr_solutions=128)
 try:r=E.enumerate(0,k,min(float(np.diag(J).min())*2.,det*.9),0)
 except EnumerationError:return math.sqrt(targets[k]/det),None
 W=np.array([np.rint(rr[1]).astype(np.int64)for rr in r]);W=np.unique(np.vstack([np.zeros((1,k),dtype=np.int64),W,-W]),axis=0);zero=int(np.flatnonzero(np.all(W==0,axis=1))[0]);A=W@J@W.T;nn=np.diag(A);D=nn[:,None]+nn[None,:]-2*A;assert np.all(D>=0)
 best=math.sqrt(targets[k]/det);bestS=[zero];bestP=0
 # Each returned cluster has its true exact diameter; greedy choices affect only discovery.
 seeds=np.argsort(nn)[:min(36,len(W))]
 for seed in seeds:
  idsW=[zero] if seed==zero else[zero,int(seed)];used=np.zeros(len(W),dtype=bool);used[idsW]=True;cost=D[idsW].max(axis=0);diam=int(cost[idsW].max())
  for size in range(len(idsW),min(18,len(W))+1):
   if size>len(idsW):
    choose=cost.copy();choose[used]=np.iinfo(np.int64).max;t=int(np.argmin(choose));idsW.append(t);used[t]=True;diam=max(diam,int(cost[t]));cost=np.maximum(cost,D[t])
   beta=diam/det
   if beta>=6:break
   val=size*(1-beta/6)**(dim/2)*math.sqrt(targets[k]/det)
   if val>best+1e-10:best=val;bestS=idsW.copy();bestP=diam
 evaluations+=1
 return best,{'B':B,'K':K,'det':det,'J':J,'W':W[bestS],'diameter_numerator':bestP}

def certify(data,ratio,k):
 B=data['B'];K=data['K'];det=int(sp.Matrix(K).det());assert det==data['det'];S=smith_normal_form(sp.Matrix(B),domain=sp.ZZ);idx=math.prod(abs(int(S[i,i]))for i in range(k))
 if idx!=1:return False
 W=data['W'];J=data['J'];assert np.array_equal(K@J,det*np.eye(k,dtype=np.int64));D=np.array([int((x-y)@J@(x-y))for x in W for y in W]);beta=Fraction(int(D.max()),det);m=len(W);dim=48-k;rat2=Fraction(m*m,det)*((6-beta)/4)**dim/target2[k]
 assert rat2>1
 row={'dimension':dim,'normal_rank':k,'normal_determinant':det,'saturation_index':idx,'cluster_size':m,'cluster_diameter_squared':str(beta),'density_ratio_squared':str(rat2),'density_ratio':math.sqrt(float(rat2)),'center_density':math.sqrt(float(target2[k]*rat2)),'seconds':time.time()-start,'minimum_dependency':'published minimum6 ofP48n, not re-enumerated here'}
 np.savez(P/f'{a.tag}_candidate_dim{dim}.npz',G=G,B=B,K=K,dual_integer_coordinates=W,adjugate=J);(P/f'{a.tag}_candidate_dim{dim}.json').write_text(json.dumps(row,indent=2));history.append(row);print('CANDIDATE',json.dumps(row),flush=True);return True

seeds={}
for k in wanted:
 seeds[k]=[]
 for mat in seedrows[k]:
  B=mat.astype(np.int16).copy();B*=np.sign(B[np.arange(k),np.argmax(B!=0,axis=1)])[:,None];seeds[k].append([lookup[v.tobytes()]for v in B])
for k in wanted:
 val,data=evaluate(seeds[k][0]);print('CONTROL',k,val,'cluster',None if data is None else len(data['W']),flush=True)
 if k==5:assert val>0.999999,'Known43D cluster not recovered by discovery control'
step=0;chains={k:(seeds[k][0],1.)for k in targets};bestweak={k:0. for k in wanted}
while time.time()-start<a.seconds:
 k=wanted[step%len(wanted)];ids,old=chains[k];cand=ids.copy();remove=int(rng.choice([1,1,1,2,2,3]));drop=set(map(int,rng.choice(k,remove,replace=False)));cand=[v for i,v in enumerate(cand)if i not in drop]
 for rep in range(remove):
  q=scores(cand);qmin=q.min()
  if step%5==0:
   options=np.flatnonzero(q<=qmin*(1.4 if step%20 else 2.0)+1e-7)
  else:options=np.flatnonzero(q<qmin+1e-7)
  cand.append(int(rng.choice(options)))
 val,data=evaluate(cand)
 if val>bestweak[k]:
  bestweak[k]=val;print('SCREEN',step,k,val,'det',None if data is None else data['det'],'cluster',None if data is None else len(data['W']),'seconds',time.time()-start,flush=True)
 if val>record[k]+1e-8 and data is not None:
  if certify(data,val,k):record[k]=val
 temperature=.035*(.3+.7*((step%160)/160))
 if val>=old or rng.random()<math.exp(max(-100,(math.log(max(val,1e-100))-math.log(old))/temperature)):chains[k]=(cand,val)
 if step%160==159:chains={kk:(seeds[kk][int(rng.integers(len(seeds[kk])))].copy(),1.)for kk in targets}
 step+=1
 if step%40==0:print('STEP',step,'chains',{kk:round(x[1],5)for kk,x in chains.items()},'records',record,'seconds',time.time()-start,flush=True)
(P/f'{a.tag}_summary.json').write_text(json.dumps({'steps':step,'evaluations':evaluations,'best_ratios':record,'best_screened_ratios':bestweak,'seconds':time.time()-start,'new_verified_candidates':history,'complete':False},indent=2));print('DONE',step,record,time.time()-start,flush=True)
