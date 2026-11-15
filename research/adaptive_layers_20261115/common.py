"""Exact small metric audits and shared pinned-baseline arithmetic."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
from itertools import product
from math import isqrt,lcm
import json,sys
import numpy as np
from flint import fmpq_mat
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cross_dimension_lower_20261115'))
from verify_counting import theta,root_gram,determinant,decsqrt


def checkpoint():return {a['dimension']:a for a in json.loads((P/'checkpoint.json').read_text())['rows']}

def gain(ratio2):
    with localcontext()as ctx:
        ctx.prec=65
        value=100*((Decimal(ratio2.numerator)/ratio2.denominator).sqrt()-1)
        k=int(value*10**6)
    lo,hi=F(k,10**6),F(k+1,10**6)
    assert (1+lo/100)**2<=ratio2<(1+hi/100)**2
    return str(value),[str(lo),str(hi)]

def finish(n,name,d2,kind,extra=None):
    before=F(checkpoint()[n]['density_squared']);ratio2=d2/before
    assert ratio2>1
    value,interval=gain(ratio2)
    out={'dimension':n,'name':name,'density_squared':str(d2),'center_density':decsqrt(d2),
         'kind':kind,'previous_density_squared':str(before),'ratio_squared':str(ratio2),
         'gain_over_checkpoint_percent':value,'gain_over_checkpoint_interval':interval}
    if extra:out.update(extra)
    return out

def metric_audit(gram,targets):
    m=fmpq_mat([[str(F(x))for x in v]for v in gram]);r=m.nrows()
    assert len(targets)==2**r
    assert all(fmpq_mat([[str(m[i,j])for j in range(k)]for i in range(k)]).det()>0 for k in range(1,r+1))
    inv=m.inv();maximum=max(map(F,targets));bounds=[]
    for i in range(r):
        x=maximum*F(str(inv[i,i]));b=isqrt(x.numerator//x.denominator)+1
        assert b*b>x;bounds.append(b)
    den=lcm(*(F(str(m[i,j])).denominator for i in range(r)for j in range(r)))
    num=np.array([[int(F(str(m[i,j]))*den)for j in range(r)]for i in range(r)],dtype=np.int64)
    assert r*r*int(abs(num).max())*max(bounds)**2<2**63
    v=np.array(list(product(*(range(-b,b+1)for b in bounds))),dtype=np.int64)
    norms=np.einsum('vi,ij,vj->v',v,num,v)
    bits=np.sum((v%2)*(1<<np.arange(r,dtype=np.int64)),axis=1)
    nonzero=np.any(v,axis=1);minima={}
    for b,t in enumerate(targets):
        mask=nonzero&(bits==b);assert mask.any()
        q=F(int(norms[mask].min()),den);assert q>=F(t),(b,q,t)
        minima[str(b)]=str(q)
    return {'gram':[[str(F(str(m[i,j])))for j in range(r)]for i in range(r)],
            'determinant':str(m.det()),'enumeration_coefficient_bounds':bounds,
            'vectors_checked':int(nonzero.sum()),'minimum_by_parity':minima,
            'all_infinite_layer_inequalities_certified':True}
