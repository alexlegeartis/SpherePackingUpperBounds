"""Replay valid but noncompetitive16D/22D examples, including exact Schur.

These examples demonstrate applicability; they do not improve the cited bounds.
Assertions must remain enabled. Requires the package's verification dependencies.
"""
import json,subprocess,sys
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_mat as Mat,arb,ctx
assert __debug__
ctx.prec=256
root=Path(__file__).resolve().parent
cases=[(16,'candidate16_wide_guarded.json','verify_guarded.py','0.025150156','0.02492121'),(22,'witness22_degree39/dim22.json','verify.py','0.003908297','0.0037259419689206')]
for n,filename,verifier,upper,comparison in cases:
    path=root/filename
    subprocess.run([sys.executable,str(root/verifier),str(path)],check=True)
    r=json.loads(path.read_text());A=Mat([[Q(z)for z in row]for row in r['matrices'][0]]);d=A.nrows()-1
    B=Mat([[A[i+1,j+1]for j in range(d)]for i in range(d)]);b=Mat([[A[i+1,0]]for i in range(d)])
    v=B.inv()*b;assert B*v==b
    alpha=A[0,0]-(b.transpose()*v)[0,0];assert alpha>0
    # The replacement A00=b^T B^-1 b has exact PSD factorization
    # [v^T;I] B [v,I]. No pair, triple, or exterior expression changes.
    V=Mat([[v[j,0]for j in range(d)]]+[[Q(int(i==j))for j in range(d)]for i in range(d)])
    saturated=V*B*V.transpose()
    assert saturated[0,0]==A[0,0]-alpha
    assert all(saturated[i,j]==A[i,j]for i in range(d+1)for j in range(d+1)if i+j>0)
    N=Q(r.get('reciprocal_gaussian_guards',{}).get('corrected_normalization',r['normalization']))
    objective=Q(r['exact_objective'])-alpha/N
    a=arb(n)/2;s=Q(str(Fraction(r['exclusion_squared'])))
    density=arb.pi()**a/(a+1).gamma()*arb(s/4)**a*arb(objective)
    assert density<arb(upper) and density>arb(comparison)
    print('Exact Schur factorization PASS',n,'density',density,'strict upper',upper,flush=True)
print('PASS: valid but NONCOMPETITIVE16D and22D examples; no record claim.')
