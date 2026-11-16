"""Reproduce the finite simplex-normal discovery search; no minimum proof.

Use --output NEW_DIR. The exact packing audit is verify.py, independently
of this floating discovery search. Complete defaults take about two minutes.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,time
import numpy as np
from flint import fmpz_mat,fmpq_mat

P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--trials',type=int);args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    G=fmpz_mat(np.load(P.parent/'mw92_sections_20261115/parent_certificate.npz',allow_pickle=False)['G'].tolist())
    Hq=47*G.inv();assert all(x.denominator==1 for row in Hq.tolist()for x in row)
    H=np.array([[int(x)for x in row]for row in Hq.tolist()],dtype=np.int64)
    maps=[]
    for raw in np.load(P.parent/'cm58_20261116/frontier/isometries92.npz',allow_pickle=False)['maps']:
        C=fmpz_mat(raw.tolist());assert C*G*C.transpose()==G and abs(C.det())==1
        A=C.inv().transpose();assert all(x.denominator==1 for row in A.tolist()for x in row)
        maps.append(np.array([[int(x)for x in row]for row in A.tolist()],dtype=np.int64))
    def key(v):return tuple(map(int,v if v[np.flatnonzero(v)[0]]>0 else-v))
    seed=np.array(json.loads((P/'discovery_seed.json').read_text())['label'],dtype=np.int64)
    assert int(seed@H@seed)==244
    vectors=[np.array(key(seed),dtype=np.int64)];seen={key(seed)};pos=0
    while pos<len(vectors):
        v=vectors[pos];pos+=1
        for A in maps:
            w=v@A;t=key(w)
            if t not in seen:seen.add(t);vectors.append(np.array(t,dtype=np.int64))
        assert len(vectors)<=200000
    V=np.array(vectors);assert len(V)==155664
    assert 92**2*int(abs(V).max())**2*int(abs(H).max())<2**52
    assert np.all(np.einsum('vi,ij,vj->v',V,H,V)==244)
    ips=V@H@V[0];rng=np.random.default_rng(202611168688);inverse=G.inv()
    old={a['dimension']:F(a['best_checked']['density_squared'])
         for a in json.loads((P/'checkpoint.json').read_text())['rows']}
    parent=F(2**184,47**14);best={};start=time.monotonic()
    for k in [3,4,5,6]:
        n=92-k;target=-244/k;band=18 if k==3 else 22
        ids=np.flatnonzero((abs(ips-target)<=band)|(abs(-ips-target)<=band));ids=ids[ids!=0]
        if len(ids)>6500:ids=rng.choice(ids,6500,replace=False)
        W=np.vstack([V[0],V[ids]*np.where(ips[ids]>0,-1,1)[:,None]])
        Z=np.rint((W.astype(float)@H.astype(float))@W.astype(float).T).astype(np.int16)
        assert np.all(np.diag(Z)==244)
        edges=abs(Z-target)<=band;np.fill_diagonal(edges,False)
        neighbors=[np.flatnonzero(row).astype(np.int32)for row in edges]
        threshold=math.log(float(old[n]/parent));trials=args.trials or (60000 if k<=5 else 100000)
        for trial in range(trials):
            chosen=[0];options=neighbors[0]
            while len(chosen)<k and len(options):
                if len(chosen)<2:j=int(rng.choice(options))
                else:
                    sample=rng.choice(options,min(5,len(options)),replace=False)
                    errors=np.sum((Z[np.ix_(sample,chosen)]-target)**2,axis=1)
                    j=int(sample[np.argmin(errors)])
                chosen.append(j);options=np.intersect1d(options,neighbors[j],assume_unique=True)
            if len(chosen)<k:continue
            K=Z[np.ix_(chosen,chosen)].astype(float)/47;sign,logdet=np.linalg.slogdet(K)
            if sign<=0:continue
            inv=np.linalg.inv(K);beta=max(np.diag(inv).max(),
                max(inv[i,i]+inv[j,j]-2*inv[i,j]for i in range(k)for j in range(i)))
            if beta>=16:continue
            score=2*math.log(k+1)-k*math.log(4)-logdet+n*math.log1p(-beta/16)
            if score<=threshold+1e-10:continue
            B=fmpz_mat(W[chosen].tolist());A,_=B.transpose().hnf(transform=True)
            if abs(fmpz_mat(A.tolist()[:k]).det())!=1:continue
            K=fmpq_mat(B)*inverse*fmpq_mat(B.transpose());inv=K.inv()
            beta=max(F(str(inv[i,i]))for i in range(k))
            for i in range(k):
                for j in range(i):beta=max(beta,F(str(inv[i,i]+inv[j,j]-2*inv[i,j])))
            density=(k+1)**2*((16-beta)/4)**n/F(str(G.det()*K.det()))
            if density<=old[n]or(n in best and density<=F(best[n]['density_squared'])):continue
            threshold=score;best[n]={'dimension':n,'normal_labels':B.tolist(),
                'cluster_labels':[[0]*k]+[[int(i==j)for j in range(k)]for i in range(k)],
                'density_squared':str(density)}
            print(n,'trial',trial,'center density',float(density)**.5,flush=True)
            (out/'discovery.json').write_text(json.dumps({'rows':[best[x]for x in sorted(best)]},indent=2,default=int)+'\n')
        print(n,'complete finite search; seconds',time.monotonic()-start,flush=True)


if __name__=='__main__':main()
