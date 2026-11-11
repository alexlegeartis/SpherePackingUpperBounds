"""Independent exact angular signs and full cyclic domain in (u,v=2ut).

Unlike the producer's (t,u) rectangle cover, this uses a rectangular domain
1<=u<=5, 0<=v<=1, with no omitted sorted-triangle rectangles.
"""
import hashlib,json,math,sys
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as Poly, arb,arb_poly,ctx
if not __debug__: raise RuntimeError('Assertions must be enabled')
ctx.prec=320
D=Path(__file__).resolve().parent
n=int(sys.argv[1]); assert n in (22,23)
hashes=json.loads((D/'source_hashes.json').read_text())
file=f'angular{n}.json'
assert hashlib.sha256((D/file).read_bytes()).hexdigest()==hashes[file]['sha256']
d=json.loads((D/file).read_text())
assert d['dimension']==n and d['degree']==14
p=Poly(list(map(Q,d['polynomial_coefficients']))); deg=p.degree(); assert deg==14
B=Q(d['B_exact']); c=Q(d['cutoff']); assert c==Q(5001,10000) and p(1)==B
def chain(p):
    seq=[p,p.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]); assert rem
        z=rem.numer(); seq.append(Poly(z/abs(z.content())))
    return seq
def variation(seq,t):
    vals=[p(t) for p in seq]; signs=[int(x>0)-int(x<0) for x in vals if x]
    return sum(x!=y for x,y in zip(signs,signs[1:]))
seq=chain(p)
vleft=variation(seq,Q(-1)); vcut=variation(seq,c); vright=variation(seq,Q(1))
assert p(-1)<0 and p(c)<0 and p(1)>0 and vleft==vcut and vcut-vright==1
cap=p+Q(1,100000); seq=chain(cap)
assert cap(-1)<0 and cap(Q(1,2))<0 and variation(seq,Q(-1))==variation(seq,Q(1,2))

# First-term maximum from exact Bernstein coefficients in the actual degree14.
cache={}
def first_bound(lo,hi):
    if (lo,hi) not in cache:
        powers=[sum((p[k]*Q(math.comb(k,j))*lo**(k-j)*(hi-lo)**j for k in range(j,deg+1)),Q(0)) for j in range(deg+1)]
        cache[lo,hi]=max(sum((powers[j]*Q(math.comb(i,j),math.comb(deg,j)) for j in range(i+1)),Q(0)) for i in range(deg+1))
    return cache[lo,hi]
derivatives=[]; z=p
for j in range(deg+1):
    derivatives.append(arb_poly([arb(x) for x in (z/math.factorial(j)).coeffs()]))
    z=z.derivative()
def angular_enclosure(interval):
    # Actual cosine is in [-1,1], so intersect with that geometrical range.
    lo=max(arb(-1),interval.lower()); hi=min(arb(1),interval.upper())
    assert lo<=hi
    mid=(lo+hi)/2; radius=(hi-lo)/2
    center=derivatives[0](mid)
    error=sum((abs(derivatives[j](mid))*radius**j for j in range(1,deg+1)),arb(0))
    return (center-error).union(center+error)

stack=[(Q(1),Q(5),Q(0),Q(1),0,0)]
passed=depthmax=0; digest=hashlib.sha256()
while stack:
    ul,uh,vl,vh,du,dv=stack.pop()
    U=arb(ul).union(arb(uh)); V=arb(vl).union(arb(vh))
    T=V/(2*U); W=(1+U*U-V).sqrt()
    bound=first_bound(vl/(2*uh),vh/(2*ul))
    z1=(1-V/2)/W; z2=(U-T)/W
    gap1=arb(1-vh).union(arb(1-vl))
    gap2=arb(ul*ul-vh).union(arb(uh*uh-vl))
    total=arb(bound)+angular_enclosure(z1)*(-40*gap1).exp()+angular_enclosure(z2)*(-40*gap2).exp()
    if total<0:
        passed+=1; depthmax=max(depthmax,du,dv)
        digest.update(f'{ul}|{uh}|{vl}|{vh}\n'.encode()); continue
    assert max(du,dv)<27,(n,ul,uh,vl,vh)
    if du<=dv:
        um=(ul+uh)/2
        stack.extend([(ul,um,vl,vh,du+1,dv),(um,uh,vl,vh,du+1,dv)])
    else:
        vm=(vl+vh)/2
        stack.extend([(ul,uh,vl,vm,du,dv+1),(ul,uh,vm,vh,du,dv+1)])
    if passed and passed%10000==0: print(n,'angular leaves',passed,flush=True)
assert Q(40**30,math.factorial(30))>2*B*100000
out={'pass':True,'dimension':n,'angular_degree':deg,'exact_Sturm_variations':[vleft,vcut,vright],
     'cap_margin':'1/100000','domain':'1<=u<=5, 0<=v=2ut<=1',
     'interval_precision_bits':ctx.prec,'certified_rectangles':passed,'excluded_rectangles':0,
     'max_depth':depthmax,'partition_sha256':digest.hexdigest(),
     'analytic_exponential_tail_bound':True,'single_sign_change_for_scale_lemma':True}
(D/f'angular{n}_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(out,flush=True)
