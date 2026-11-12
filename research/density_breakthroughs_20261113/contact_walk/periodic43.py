#!/usr/bin/env python3
"""Full first-order periodic metric/translate cone at a published 43D record.

The 946 metric and 215 relative-translate variables are independent. Contacts
come from individually exact norm-six ambient vectors; the pool is incomplete.
A numerical descent is a lead, not a packing certificate. A numerical optimum
zero is a screen, not an exact eutaxy certificate.
"""
import json
import re
import time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import null_space
from scipy.optimize import linprog

P = Path(__file__).resolve().parent
BASE = P.parent
CERT = BASE/'seed_deformations/external_sphere_packing/packings_out/records/dim43-delta226/P48n-081e3830.cert.txt'

def mat(text, key):
    a = text.index('[[', text.index(key))
    b = text.index(']]',a)+2
    return np.array([[int(x) for x in re.findall(r'-?\d+',line)] for line in text[a:b].splitlines()],dtype=np.int64)

def main():
    start=time.time()
    text=CERT.read_text()
    G=mat(text,'ambient Gram G')
    T=mat(text,"simplex-basis coefficient rows T'")
    K=T@G@T.T
    assert sp.Matrix(K).det()==3072
    Ki=sp.Matrix(K).inv()
    Adj=np.array(Ki*3072,dtype=np.int64)
    assert np.array_equal(K@Adj,3072*np.eye(5,dtype=np.int64))
    data=np.load(BASE/'gluing/P48n_glue_pool.npz')
    assert np.array_equal(data['G'],G)
    V=data['V']
    assert np.all(np.einsum('ij,jk,ik->i',V,G,V)==6)
    S=np.vstack([np.zeros(5,dtype=np.int64),np.eye(5,dtype=np.int64)])
    normal=V@G@T.T
    L=np.linalg.cholesky(G)
    E=null_space(T@L)
    assert E.shape==(48,43)
    assert np.linalg.norm(E.T@E-np.eye(43))<1e-10
    projected=L@E
    contacts=[];pairrows=[];pairs=[]
    for i in range(6):
        for j in range(i+1,6):
            u=S[i]-S[j]
            norm_num=int(u@Adj@u)
            if norm_num!=896:
                assert norm_num<896
                continue
            yes=np.all(normal==u,axis=1)
            neg=np.all(normal==-u,axis=1)
            X=np.vstack([V[yes],-V[neg]])
            W=X@projected
            assert np.max(abs(np.einsum('ij,ij->i',W,W)-137/24))<1e-8
            contacts.extend(X)
            pairrows.extend([(i,j)]*len(X))
            pairs.append({'i':i,'j':j,'contacts':len(X),'beta':'7/24'})
    contacts=np.array(contacts,dtype=np.int64)
    pairrows=np.array(pairrows,dtype=np.int64)
    W=contacts@projected
    n=43;ii,jj=np.triu_indices(n);f=np.where(ii==jj,1.,2.)
    A=np.zeros((len(W),len(ii)+5*n))
    A[:,:len(ii)]=W[:,ii]*W[:,jj]*f
    for i in range(1,6):
        A[:,len(ii)+(i-1)*n:len(ii)+i*n]=2*W*( (pairrows[:,0]==i).astype(int)-(pairrows[:,1]==i).astype(int))[:,None]
    objective=np.zeros(A.shape[1]);objective[:len(ii)]=(ii==jj)
    print('contacts',len(W),'variables',A.shape[1],'pairs',pairs,flush=True)
    rng=np.random.default_rng(4313)
    selected=rng.choice(len(W),min(2500,len(W)),replace=False)
    history=[]
    for iteration in range(30):
        out=linprog(objective,A_ub=-A[selected],b_ub=np.zeros(len(selected)),bounds=(-1,1),method='highs-ipm',options={'threads':1,'time_limit':90})
        rec={'iteration':iteration,'success':bool(out.success),'message':out.message,'constraints':len(selected),'seconds':time.time()-start}
        if not out.success:
            history.append(rec);break
        slopes=A@out.x
        bad=np.flatnonzero(slopes<-1e-8)
        rec.update({'trace':float(out.fun),'violated_other_contacts':len(bad),'minimum_slope':float(slopes.min())})
        history.append(rec);print(json.dumps(rec),flush=True)
        if len(bad):
            selected=np.unique(np.concatenate([selected,bad[np.argsort(slopes[bad])[:1500]]]))
        else:
            np.savez_compressed(P/'periodic43_cone.npz',G=G,T=T,K=K,E=E,projected_basis=projected,contacts=contacts,pairs=pairrows,selected=selected,direction=out.x,dual=out.ineqlin.marginals)
            break
    report={'dimension':43,'seed_center_density':float(6*(137/96)**21.5/np.sqrt(3072)),'cosets':6,'metric_variables':len(ii),'translate_variables':215,'ambient_pool_count':len(V),'contact_count':len(W),'active_pairs':pairs,'history':history,'scope':'floating full periodic first-order cone; exact individual contact norms; no new packing or exact local-optimum claim'}
    (P/'periodic43_report.json').write_text(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
