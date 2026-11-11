"""Read-only exact origin-saturation analysis; does not register a certificate.
Requires python-flint only for final outward interval arithmetic.
Source full-domain audits are separate proof dependencies.
"""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path
from flint import arb, fmpq, ctx
if not __debug__: raise RuntimeError('Assertions must be enabled')
ctx.prec=384

def solve(A,b):
    M=[list(row)+[v] for row,v in zip(A,b)]
    for i in range(len(M)):
        assert M[i][i]>0
        for j in range(i+1,len(M)):
            t=M[j][i]/M[i][i]
            for k in range(i,len(M)+1): M[j][k]-=t*M[i][k]
    x=[F(0)]*len(M)
    for i in range(len(M)-1,-1,-1):
        x[i]=(M[i][-1]-sum(M[i][j]*x[j] for j in range(i+1,len(M))))/M[i][i]
    return x

def audit(path,threshold=None):
    data=path.read_bytes();j=json.loads(data);A=[[F(z) for z in row] for row in j['matrices'][0]]
    B=[row[1:] for row in A[1:]];b=A[0][1:];v=solve(B,b)
    assert all(sum(B[i][k]*v[k] for k in range(len(b)))==b[i] for i in range(len(b)))
    alpha=A[0][0]-sum(x*y for x,y in zip(b,v));assert alpha>0
    W=[v]+[[F(i==k) for k in range(len(b))] for i in range(len(b))]
    sat=[[sum(W[i][p]*B[p][q]*W[k][q] for p in range(len(b)) for q in range(len(b))) for k in range(len(A))] for i in range(len(A))]
    assert all(sat[i][k]==A[i][k]-(alpha if i==k==0 else 0) for i in range(len(A)) for k in range(len(A)))
    norm=F(j['normalization']);assert norm>0
    objective=(F(j['physical_polynomial'][0])+A[0][0])/norm;assert objective==F(j['exact_objective'])
    newobj=objective-alpha/norm
    n=j['dimension'];a=arb(n)/2;s=F(j['exclusion_squared']);vol=arb.pi()**a/(a+1).gamma()*arb(fmpq(str(s/4)))**a
    density=vol*arb(fmpq(str(newobj)));decrement=vol*arb(fmpq(str(alpha/norm)))
    if threshold is not None: assert density<arb(threshold)
    return {'source':str(path),'sha256':hashlib.sha256(data).hexdigest(),'dimension':n,'exact_factorization':'PASS with independent Fraction elimination','alpha':str(alpha),'normalization':str(norm),'saturated_exact_objective':str(newobj),'density_interval':str(density),'density_decrement_interval':str(decrement),'proposed_outward_threshold':threshold,'registration_status':'not registered by this analysis; source full-domain verification remains an external dependency'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--threshold');args=p.parse_args();print(json.dumps(audit(args.source,args.threshold),indent=2))
