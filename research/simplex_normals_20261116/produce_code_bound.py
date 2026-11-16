"""Produce an exact Delsarte LP certificate for an even binary code.

This optional rational simplex solver is discovery/reproduction machinery.
verify_craig.py checks the saved polynomial certificates without a solver.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse,json
from flint import fmpq_mat
from verify_craig import kw,polynomial_bound


def solve(N,d,shortened=0):
    n=N-shortened;weights=list(range(d,n+1,2));m=len(weights);r=n//2
    assert 0<d<=n and d%2==0 and 0<=shortened<N
    basis=list(range(m,m+r));A=fmpq_mat(r+1,m+r+1)
    for j in range(r):
        for i,w in enumerate(weights):A[j,i]=str(-F(kw(n,j+1,w),comb(n,j+1)))
        A[j,m+j]=1;A[j,m+r]=1
    for i in range(m):A[r,i]=-1
    # Maximize sum(distance-distribution entries), beginning with the
    # feasible one-word code. Bland's rule handles entering/leaving ties.
    for _ in range(10000):
        col=next((i for i in range(m+r)if A[r,i]<0),None)
        if col is None:break
        choices=[j for j in range(r)if A[j,col]>0];assert choices
        row=min(choices,key=lambda j:(A[j,m+r]/A[j,col],basis[j]))
        pivot=A[row,col];R=fmpq_mat([[A[row,i]/pivot for i in range(m+r+1)]])
        C=fmpq_mat([[A[j,col]]for j in range(r+1)]);A=A-C*R
        for i in range(m+r+1):A[row,i]=R[0,i]
        basis[row]=col
    else:raise RuntimeError('Exact simplex pivot cap reached.')
    y=[F(str(A[r,m+j]))for j in range(r)];assert all(t>=0 for t in y)
    degrees=[j+1 for j,t in enumerate(y)if t]
    result={'shortened':shortened,'length_after_shortening':n,
            'krawtchouk_degrees':[0]+degrees,
            'normalized_coefficients':['1']+[str(y[j-1])for j in degrees],
            'exact_code_size_bound':str((1+sum(y))*2**shortened)}
    polynomial_bound(result,N,d)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--length',type=int,required=True)
    ap.add_argument('--distance',type=int,required=True);ap.add_argument('--shortened',type=int,default=0)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    result=solve(args.length,args.distance,args.shortened)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('Exact even-code size bound',result['exact_code_size_bound'])
