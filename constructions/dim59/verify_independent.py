"""Second exact implementation: FLINT HNF kernels and Sylvester determinants.

This does not import generate.py or verify.py. It also reconstructs the 59D
Gamma72 comparison period, retaining the two published minimum theorems.
"""
from pathlib import Path
import hashlib
import html
import itertools
import json
import re
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpq

if not __debug__:
    raise RuntimeError("Exact proof checks require Python assertions; run without -O.")

P = Path(__file__).resolve().parent
ROOT = P.parents[1]


def integer(rows): return fmpz_mat([[int(x) for x in row] for row in rows])
def rational(M): return fmpq_mat(M.tolist())
def eye(n): return fmpz_mat([[int(i==j) for j in range(n)]for i in range(n)])
def read(name):
    return integer([[int(v)for v in line.split()]for line in(P/"data"/name).read_text().splitlines()])
def rows(M,a,b): return type(M)(M.tolist()[a:b])


def positive(M):
    values = M.tolist()
    for k in range(1,M.nrows()+1):
        assert integer([row[:k]for row in values[:k]]).det()>0


def kernel(A):
    H,T = A.transpose().hnf(transform=True)
    k,n = A.nrows(),A.ncols()
    assert T*A.transpose()==H and abs(T.det())==1
    assert H==integer([[int(i==j)for j in range(k)]for i in range(n)])
    Ti = integer(T.inv().tolist())
    assert T*Ti==eye(n)
    return rows(T,k,n),rows(T,0,k),Ti


def main():
    source=ROOT/"research/novel_routes_20261113/dual64/CQ64b.html"
    raw=source.read_text().split('<a NAME="GRAM"',1)[1].split('<a NAME="DIVISORS"',1)[0]
    raw=html.unescape(re.sub(r'<[^>]*>',' ',raw))
    nums=list(map(int,re.findall(r'-?\d+',raw)))
    assert nums[:2]==[64,0]and len(nums)==2082
    values=[[0]*64 for _ in range(64)];pos=2
    for i in range(64):
        for j in range(i+1):values[i][j]=values[j][i]=nums[pos];pos+=1
    G=integer(values);assert G.det()==3**32
    positive(G)
    Gi=G.inv();A=read("normal_coefficients.txt");S=read("dual_cluster.txt")
    assert A.nrows()==5 and A.ncols()==64
    assert S==integer([[0]*5]+eye(5).tolist())
    K=A*Gi*A.transpose();assert 3*K==read("normal_gram_numerator.txt")
    assert K.det()==fmpq(33800,81)
    N2,right,Ti=kernel(A)
    Q2=N2*G*N2.transpose()
    assert Q2.det()==33800*3**28
    N1=read("section_basis.txt");Q1=read("period_gram.txt")
    change=N1*Ti
    assert all(change[i,j]==0 for i in range(59)for j in range(5))
    W=integer([row[5:]for row in change.tolist()])
    assert abs(W.det())==1 and W*N2==N1
    assert Q1==W*Q2*W.transpose()
    lifts1=read("ambient_lifts.txt");lifts2=S*right
    assert lifts1*A.transpose()==lifts2*A.transpose()==S
    C=rational(read("shift_numerators.txt"))/int(read("shift_denominator.txt")[0,0])
    assert C*N1==rational(lifts1)-S*K.inv()*A*Gi
    distances=[]
    for i,j in itertools.combinations(range(6),2):
        d=integer([[S[i,c]-S[j,c]for c in range(5)]])
        distances.append((d*K.inv()*d.transpose())[0,0])
    beta=max(distances);assert beta==fmpq(45,104)
    minimum=12-beta;assert minimum==fmpq(1203,104)
    contact=json.loads((P/"data/contact.json").read_text())
    x=integer([contact["ambient_coefficients"]]);i,j=contact["i"],contact["j"]
    d=integer([[S[i,c]-S[j,c]for c in range(5)]])
    assert x*A.transpose()==d and(x*G*x.transpose())[0,0]==12
    projected=rational(x)-d*K.inv()*A*Gi
    assert(projected*G*projected.transpose())[0,0]==minimum
    delta2=36*(minimum/4)**59/Q2.det()
    assert delta2==fmpq(3**33*401**59,2**296*5**2*13**61)

    # Independent reconstruction of the stated comparison, not a table lookup.
    source72=ROOT/"research/novel_routes_20261113/higher_parent/neb72.html"
    raw=source72.read_text().split('<a NAME="GRAM"',1)[1].split('<a NAME="DET"',1)[0]
    nums=list(map(int,re.findall(r'-?\d+',re.sub(r'<[^>]*>',' ',raw))))
    assert nums[0]==72 and len(nums)==5185
    G72=integer([nums[1+72*i:1+72*(i+1)]for i in range(72)])
    assert G72==G72.transpose() and G72.det()==1
    positive(G72)
    baseline=ROOT/"research/novel_routes_20261113/audits/gamma72_leech_section_rank13.npz"
    b=np.load(baseline);B=integer(b["B"].tolist())
    assert G72==integer(b["G"].tolist())
    K72=B*G72*B.transpose();assert K72==integer(b["K"].tolist())
    assert K72.det()==2**23
    baseN,_,_=kernel(B*G72)
    assert baseN*G72*B.transpose()==integer([[0]*13 for _ in range(59)])
    baseQ=baseN*G72*baseN.transpose();assert baseQ.det()==2**23
    baseline2=fmpq(2**59,baseQ.det());assert baseline2==2**36
    assert delta2>baseline2
    report={"all_checks_passed":True,"implementation":"python-flint; independent HNF kernel and Sylvester positivity",
        "dimension":59,"cosets":6,"minimum_squared_exact":str(minimum),
        "center_density_squared_exact":str(delta2),"period_basis_change_determinant":int(W.det()),
        "comparison_reconstructed_exactly":True,"comparison_center_density_squared":str(baseline2),
        "comparison_kind":"internal_repository_comparison",
        "comparison_is_published_59D_record":False,
        "inherited_minimum_theorems":{"Ne64":12,"Gamma72":8},
        "literature_status_report":"../../research/covering_and_deep_cosets_20261115/PRIORITY.md",
        "data_sha256":{f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in sorted((P/"data").iterdir())if f.is_file()},
        "source_sha256":{"Ne64":hashlib.sha256(source.read_bytes()).hexdigest(),"Gamma72":hashlib.sha256(source72.read_bytes()).hexdigest()},
        "comparison_certificate_sha256":hashlib.sha256(baseline.read_bytes()).hexdigest()}
    (P/"independent_audit.json").write_text(json.dumps(report,indent=2)+"\n")
    print("PASS: independent exact 59D packing and Gamma72 comparison; inherited minima 12 and 8.")


if __name__=="__main__":main()
