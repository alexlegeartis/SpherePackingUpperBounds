"""Exact upper witnesses for failed holes; no inference from unresolved rows."""
from pathlib import Path
import json
import numpy as np
from flint import fmpz_mat

if not __debug__:
    raise RuntimeError("Run the exact witness audit without Python -O.")

P=Path(__file__).resolve().parent
Q=np.load(P.parent/"novel_routes_20261113/higher_parent/ne64_primal_best.npz")["R"].astype(np.int64)
reports=[]
for name,modulus,cut in (("half_coset_rejections.npz",2,34),("third_coset_screen.npz",3,66)):
    a=np.load(P/name)
    C=a["residues"].astype(np.int64)
    norms=a["norms"].astype(np.int64)
    known=norms>=0
    Y=a["witnesses"][known].astype(np.int64)
    assert np.all((Y-C[known])%modulus==0)
    assert 64**2*int(abs(Q).max())*int(abs(Y).max())**2<2**63
    assert np.array_equal(np.einsum('vi,ij,vj->v',Y,Q,Y),norms[known])
    assert np.all(norms[known]<=cut)
    if modulus==2:
        assert len(C)==100 and np.all(known)
    else:
        assert len(C)==10000
        assert np.all((C@Q)%3==0)
        assert np.all(np.einsum('vi,ij,vj->v',C,Q,C)%18==12)
    reports.append({"file":name,"classes":len(C),"exact_rejections":int(known.sum()),
                    "unresolved":int((~known).sum()),"inclusive_numerator_cutoff":cut,
                    "modulus":modulus})
abstract=json.loads((P/"gamma72_66_abstract.json").read_text())
K=fmpz_mat(abstract["normal_gram"]);H=K.inv();S=fmpz_mat(abstract["dual_labels"])
assert K.det()==21875
for j in range(1,7):assert fmpz_mat([row[:j]for row in K.tolist()[:j]]).det()>0
diameter=0
for i in range(S.nrows()):
    for j in range(i):
        d=fmpz_mat([[S[i,c]-S[j,c]for c in range(6)]])
        diameter=max(diameter,(d*H*d.transpose())[0,0])
assert str(diameter)=="2/5"
assert abstract["ambient_embedding_found"] is False
report={"exact_upper_witnesses_pass":True,"screens":reports,
        "abstract_66D_gram_and_diameter_pass":True,"abstract_66D_is_a_packing_certificate":False,
        "scope":"Only explicit witnesses reject classes. Unresolved/pruned/capped searches prove no hole; the unembedded normal is not a packing."}
(P/"screen_audit.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
