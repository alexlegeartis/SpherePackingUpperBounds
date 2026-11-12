#!/usr/bin/env python3
"""Generate exact minimum-vector orbits in the catalogue's 56D Delta2 form.

Norm tests and group invariance use integer arithmetic. Orbit completeness is
reported against the catalogue's count, not claimed as a new SVP proof.
"""
from pathlib import Path
import json,re,time
import numpy as np
from scipy.linalg import eigvalsh

P=Path(__file__).resolve().parent

def section(s,key):
    m=re.search(r'<a\s+NAME="'+key+r'[^\"]*">.*?</a>(.*?)(?=<p><li>|<a\s+NAME=|<hr>|$)',s,re.S|re.I)
    return re.sub(r'<[^>]+>',' ',m.group(1))

def main():
    start=time.time();s=(P/'sources/L_56a.html').read_text();n=56
    raw=list(map(int,re.findall(r'-?\d+',section(s,'GRAM'))));assert raw[:2]==[n,0]
    G=np.zeros((n,n),dtype=np.int64);i,j=np.tril_indices(n);G[i,j]=raw[2:];G[j,i]=raw[2:]
    raw=list(map(int,re.findall(r'-?\d+',section(s,'SUBGROUP_GENERATORS'))));ng=raw[0];raw=raw[1:];As=[]
    for a in range(ng):
        assert raw[:2]==[n,n];A=np.array(raw[2:2+n*n],dtype=np.int64).reshape(n,n);raw=raw[2+n*n:]
        if not np.array_equal(A@G@A.T,G):
            assert np.array_equal(A.T@G@A,G);A=A.T
        assert np.array_equal(A@G@A.T,G);As.append(A)
    assert not raw
    seen=set();vectors=[]
    def put(V):
        added=[]
        for v in V:
            nz=np.flatnonzero(v);assert len(nz)
            if v[nz[0]]<0:v=-v
            assert np.max(abs(v))<32768
            key=v.astype(np.int16).tobytes()
            if key not in seen:seen.add(key);vectors.append(v.copy());added.append(v)
        return np.array(added,dtype=np.int64).reshape(-1,n)
    current=put(np.eye(n,dtype=np.int64)[np.diag(G)==28])
    while len(current):
        current=put(np.vstack([current@A for A in As]))
        print('new',len(current),'total',len(vectors),'seconds',time.time()-start,flush=True)
        assert len(vectors)<=53940
    V=np.array(vectors,dtype=np.int64)
    assert np.all(np.einsum('bi,ij,bj->b',V,G,V)==28)
    M=V.T@V;F=G@M@G
    eig=eigvalsh(F,G)
    np.savez_compressed(P/'L56a_orbit.npz',G=G,V=V,generators=np.array(As),frame=F)
    report={'dimension':56,'vectors_up_to_sign':len(V),'catalogue_signed_kissing':107880,'frame_eigenvalues':eig.tolist(),'seconds':time.time()-start,'scope':'exact isometries and contact witnesses; no new minimum certificate'}
    (P/'L56a_orbit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

if __name__=='__main__':main()
