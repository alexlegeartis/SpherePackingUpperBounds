"""Exact coefficient and atom comparison. No sampled point decides equality."""
from fractions import Fraction as Q
import symalg as sa
from symalg import Sym


def require_equal(a, b, label):
    if not (a-b).is_zero():
        raise AssertionError(label + ': symbolic expressions differ')


def homogenize(p, degree=18):
    """(2rs)^degree P(r,s,-1+z*((r+s)^2-49/16)/(2rs)).

    The third polynomial generator is named t internally but stands for z in
    this return value. All operations are exact in Q[r,s,z].
    """
    r,s,z=sa.CTX.gens()
    denom=2*r*s
    numer=z*((r+s)**2-sa.fq(Q(49,16)))-denom
    coeffs=[{} for _ in range(degree+1)]
    for (i,j,k),value in p.to_dict().items():
        assert k <= degree, 'angular degree exceeds homogenizing degree'
        coeffs[k][(i,j,0)]=value
    return sum((sa.CTX.from_dict(row)*numer**k*denom**(degree-k)
                for k,row in enumerate(coeffs) if row), sa.ZERO_P)


def branch_poly(branch, swapped=False):
    assert type(branch['dd']) is int and branch['dd'] > 0
    assert all(len(exp)==3 and all(type(e) is int and e>=0 for e in exp)
               and type(c) is int for exp,c in branch['C'].items())
    d={((j,i,k) if swapped else (i,j,k)):sa.fq(Q(v,branch['dd']))
       for (i,j,k),v in branch['C'].items()}
    return sa.CTX.from_dict(d)


def compare(s1,s5):
    assert s1['breaks']==s5['breaks']
    assert list(map(Q,s1['breaks']))==list(map(Q,[1,'501473/500000','1005643/1000000',
         '50379/50000','50767/50000','102597/100000','3/2','7/4']))
    assert Q(s5['a'])==Q(133,250)
    assert 2+2*Q(s5['a'])>Q(49,16)
    assert s1['angular_degree']==18 and s1['pair_pieces']==7
    cells={(i,j) for i in range(5) for j in range(5)}
    assert set(s1['corrections'])==set(s5['cells'])==cells
    assert set(s1['endpoints'])==set(s5['endpoints'])==cells
    assert set(s1['branches'])=={(i,j) for i in range(5) for j in range(i,5)}
    rows=[]
    for i,j in sorted(cells):
        row=s5['cells'][i,j]
        kc,kf=Sym.fromkey(row['kcert']),Sym.fromkey(row['kformula'])
        require_equal(kc,kf,f'kcert exact-kink expression versus kformula cell {(i,j)}')
        correction=Sym.fromkey(s1['corrections'][i,j])
        positive_kernel=correction-kf
        assert positive_kernel.is_free(), f'nonpolynomial terms fail to cancel at {(i,j)}'
        p=positive_kernel.free()
        assert p.is_poly(), f'unexpected denominator at {(i,j)}'
        actual=homogenize(p.n)
        expected=branch_poly(s1['branches'][min(i,j),max(i,j)],i>j)
        assert actual==expected, f'polynomial coefficient mismatch at {(i,j)}'
        endpoint=s5['endpoints'][i,j]
        e5=Sym.fromkey(endpoint['exact'])
        require_equal(e5,Sym.fromkey(endpoint['lower'])+Sym.fromkey(endpoint['dropped']),
                      f'endpoint lower envelope cell {(i,j)}')
        ep=Sym.fromkey(s1['endpoints'][i,j])-e5
        assert ep.is_free() and ep.free().is_poly(), f'endpoint atoms differ at {(i,j)}'
        expected0=sa.CTX.from_dict({exp:c for exp,c in expected.to_dict().items() if exp[2]==0})
        assert ep.free().n*(2*sa.GR*sa.GS)**18==expected0, f'endpoint polynomial mismatch at {(i,j)}'
        # Every lower()-selected cap denominator must be c*(1-t^2), c>0;
        # this is positive on the entire open angular interval, not just a probe.
        for key in s5['lower_log'][i,j]:
            v=Sym.fromkey(key)
            assert v.is_free() and v.free().is_poly()
            d=v.free().n.to_dict(); c=d.get((0,0,0),0)
            assert c>0 and d=={(0,0,0):c,(0,0,2):-c}
        rows.append(dict(cell=[i,j],polynomial_coefficients=len(actual.to_dict()),
                         nonpolynomial_monomials=len(kf.terms)-1,
                         cap_denominator_checks=len(s5['lower_log'][i,j])))
    return rows


def negative_controls(s1,s5):
    """Corrupt fresh in-memory results, never source certificates."""
    import copy
    cases=[]
    a=copy.deepcopy(s5)
    a['cells'][0,0]['kcert']=(Sym.fromkey(a['cells'][0,0]['kcert'])+1).key()
    cases.append(('changed_stage5_constant',s1,a))
    a=copy.deepcopy(s1)
    first=next(iter(a['branches'][0,0]['C']))
    a['branches'][0,0]['C'][first]+=1
    cases.append(('changed_stage1_coefficient',a,s5))
    a=copy.deepcopy(s5)
    a['cells'].pop((4,4))
    cases.append(('missing_radial_cell',s1,a))
    a=copy.deepcopy(s5)
    a['endpoints'][0,0]['exact']=(Sym.fromkey(a['endpoints'][0,0]['exact'])+1).key()
    cases.append(('changed_endpoint_value',s1,a))
    a=copy.deepcopy(s1)
    a['branches'][0,0]['dd']=Q(3,2)
    cases.append(('nonintegral_denominator',a,s5))
    result=[]
    for name,left,right in cases:
        try:
            compare(left,right)
        except (AssertionError,KeyError):
            result.append(dict(name=name,status='REJECTED'))
        else:
            raise AssertionError('negative control accepted: '+name)
    return result
