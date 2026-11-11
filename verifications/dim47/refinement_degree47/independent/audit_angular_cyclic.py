"""Independent exact angular and complete cyclic audit; no optimization.

Rebuilds the dimension-specific harmonic basis by recurrence, uses rational
Sturm chains, and proves the triangle domain with a fresh interval subdivision.
The geometry radicand uses exact monotone endpoint bounds rather than the
producer's natural interval expression.
"""
from pathlib import Path
import json
import math
import sys
import time
from flint import fmpq as Q, fmpq_poly as P, arb, arb_poly, ctx

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0)
ctx.prec = 256
D = Path(__file__).resolve().parent
r = json.loads((D/'inputs/angular47.json').read_text())
n = r['dimension']
assert n == 47 and r['degree'] == 19 and r['alpha'] == 40
p = P(list(map(Q, r['polynomial_coefficients'])))
b = list(map(Q, r['gegenbauer_coefficients']))
B = Q(r['B_exact'])
cut = Q(r['cutoff'])
delta = Q(r['negative_cap_delta'])
assert cut == Q(5001,10000) and delta == Q(1,100000)
x = P([0,1])
G = [P([1]), x]
for ell in range(1,19):
    G.append(((2*ell+n-2)*x*G[-1]-ell*G[-2])/(ell+n-2))
assert all(g(1)==1 for g in G)
assert b[0]==1 and all(c>=0 for c in b)
assert sum((c*g for c,g in zip(b,G)), P([])) == p
assert len(b)==20 and p.degree()==19 and sum(b)==B and p(1)==B
moment=Q(1)
b0=Q(0)
for j in range(10):
    if j:
        moment *= Q(2*j-1, n+2*j-2)
    b0 += p[2*j]*moment
assert b0==1


def sturm(poly):
    chain=[poly,poly.derivative()]
    while chain[-1].degree()>0:
        rem=-(chain[-2]%chain[-1])
        assert rem, 'Repeated polynomial factor requires separate treatment'
        numer=rem.numer()
        chain.append(P(numer/abs(numer.content())))
    assert chain[-1]
    return chain


def variations(chain,point):
    vals=[g(point) for g in chain]
    sg=[1 if c>0 else -1 for c in vals if c]
    return sum(a!=b for a,b in zip(sg,sg[1:]))


chain=sturm(p)
vc=[variations(chain,t) for t in [Q(-1),cut,Q(1)]]
assert p(-1)<0 and p(cut)<0 and p(1)>0
assert vc[0]==vc[1] and vc[1]-vc[2]==1
cap=p+delta
cc=sturm(cap)
assert variations(cc,Q(-1))==variations(cc,Q(1,2))
assert cap(-1)<0 and cap(Q(1,2))<0
lo,hi=cut,Q(1)
for _ in range(96):
    mid=(lo+hi)/2
    val=p(mid)
    assert val
    if val<0:lo=mid
    else:hi=mid
assert p(lo)<0 and p(hi)>0
assert Q(40**30,math.factorial(30))>2*B/delta
angular={'pass':True,'dimension':n,'degree':p.degree(),'b0_exact':str(b0),
         'B_exact':str(B),'B_float':float(B),'harmonic_expansion_exact':True,
         'constant_coefficient_spherical_moment_exact':True,
         'all_harmonic_coefficients_nonnegative':True,
         'sturm_variations_minus_one_cut_one':vc,
         'strict_cap_margin_exact':str(delta),
         'unique_root_bracket_exact':[str(lo),str(hi)],
         'unique_root_bracket_float':[float(lo),float(hi)],
         'tail_exponential_majorant_exact':True}
(D/'angular_result.json').write_text(json.dumps(angular,indent=2)+'\n')
print('ANGULAR_PASS',json.dumps(angular),flush=True)

AP=arb_poly([arb(c) for c in p.coeffs()])
cache={}


def first_upper(a,b):
    key=(a,b)
    if key not in cache:
        translated=p(P([a,b-a]))
        bern=[sum((translated[j]*Q(math.comb(i,j),math.comb(19,j))
                   for j in range(i+1)),Q(0)) for i in range(20)]
        cache[key]=max(bern)
    return arb(cache[key])


def angular_interval(t):
    t=t.intersection(arb(0,1))
    coeff=AP(arb_poly([t.mid(),t.rad()]))
    err=sum((abs(c) for c in coeff.coeffs()[1:]),arb(0))
    return coeff[0]+arb(0,err.upper())


stack=[(Q(0),Q(1,2),Q(1),Q(5),0,0)]
passed=discarded=processed=depth=0
began=time.time()
while stack:
    a,b,c,d,lt,lu=stack.pop()
    processed+=1
    if 2*a*c>1:
        discarded+=1
        continue
    t=arb(a).union(arb(b))
    u=arb(c).union(arb(d))
    # 1+u^2-2tu increases in u and decreases in t on this rectangle.
    vlo=1+c*c-2*b*c
    vhi=1+d*d-2*a*d
    assert vlo>=1 and vhi>=vlo
    v=arb(vlo).union(arb(vhi)).sqrt()
    assert v.is_finite() and v>0
    gap1=1-2*t*u
    gap2=u*u-2*t*u
    # Only sorted triangles matter. Clipping to nonnegative gaps restricts
    # interval overestimation without discarding any feasible point.
    if gap1.lower()<0:
        gap1=arb(0).union(gap1.upper())
    if gap2.lower()<0:
        gap2=arb(0).union(gap2.upper())
    upper=first_upper(a,b)
    z1=(1-u*t)/v
    z2=(u-t)/v
    value=upper+angular_interval(z1)*(-40*gap1).exp()+angular_interval(z2)*(-40*gap2).exp()
    assert value.is_finite()
    if value<0:
        passed+=1
        depth=max(depth,lt,lu)
        continue
    assert max(lt,lu)<30, (a,b,c,d,str(value))
    if upper>=0 or lt<=lu:
        m=(a+b)/2
        stack.extend([(a,m,c,d,lt+1,lu),(m,b,c,d,lt+1,lu)])
    else:
        m=(c+d)/2
        stack.extend([(a,b,c,m,lt,lu+1),(a,b,m,d,lt,lu+1)])
    if processed%10000==0:
        print('CYCLIC_PROGRESS',processed,passed,discarded,len(stack),flush=True)

cyclic={'pass':True,'dimension':n,'precision_bits':ctx.prec,'alpha':40,
        'compact_domain':'0<=t<=1/2, 1<=u<=5, 2tu<=1',
        'certified_rectangles':passed,'discarded_rectangles':discarded,
        'processed_rectangles':processed,'max_axis_depth':depth,
        'changed_geometry_radicand':'exact monotone rational endpoint bounds',
        'gap_clipping':'nonnegative on the sorted-triangle feasible domain',
        'signed_secondary_terms':True,'tails_exact':True,
        'tail_first_gap_if_positive_exact':str(10*cut-1),
        'tail_second_gap_lower_exact':'24',
        'scale_extension':'single-sign-transition lemma; see AUDIT.md',
        'seconds':time.time()-began}
(D/'cyclic_result.json').write_text(json.dumps(cyclic,indent=2)+'\n')
print('CYCLIC_PASS',json.dumps(cyclic),flush=True)
