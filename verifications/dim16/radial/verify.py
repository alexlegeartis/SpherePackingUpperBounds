"""One exact epsilon=1e-6 c1 repair and a complete conditional pair proof.

The spherical average-count hypothesis is NOT established in this script.
Run without -O. All acceptance gates use exact rationals or Arb intervals.
"""
if not __debug__:raise RuntimeError('Assertions required')
import sys,json,time,math,hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
sys.set_int_max_str_digits(0);ctx.prec=384
D=Path(__file__).resolve().parent;O=D
start=time.monotonic()
def write(name,j):(O/name).write_text(json.dumps(j,indent=2)+'\n')
def q(x):return Q(str(F(str(x))))
def af(x):return arb(Q(x.numerator,x.denominator)) if isinstance(x,F) else arb(x)
source=json.loads((D/'conditional_refit_exact.json').read_text())
P=Poly([Q(x) for x in source['physical']]);H=Poly([Q(x) for x in source['fourier']])
P_unrepaired=P
raw=json.loads((D/'conditional_refit.json').read_text())
def binaryq(x):
    z=F(float(x));return Q(z.numerator,z.denominator)
L=[Poly([1]),Poly([8,-2])];XX=Poly([0,1])
for k in range(1,len(raw['radial'])-1):
    L.append(((2*k+8-2*XX)*L[-1]-(k+7)*L[-2])/(k+1))
L=[p/p(0) for p in L]
assert P==sum((binaryq(c)*p for c,p in zip(raw['radial'],L)),Poly())
assert H==sum(((-1)**k*binaryq(c)*p for k,(c,p) in enumerate(zip(raw['radial'],L))),Poly())
assert Q(source['C'])==binaryq(raw['C'])
assert Q(source['lambda_count_exact'])==binaryq(raw['lambda_count'])
N=H(0);assert N>0;epsilon=Q(1,10**6);X=Poly([0,1]);E1=1-X/4
# Normalize the source first, so the modified Fourier value at zero is exactly 1-epsilon.
P=P/N+epsilon*E1;H=H/N-epsilon*E1
C=Q(source['C'])/N;lam=Q(source['lambda_count_exact'])/N
assert H(0)==1-epsilon
newN=H(0);P/=newN;H/=newN;C/=newN;lam/=newN
assert H(0)==1
variant=dict(certified_bound=False,count_conditional=True,source='conditional_refit_exact.json',
    epsilon_exact=str(epsilon),repair='Normalize source by exact H(0), add epsilon to Laguerre c1, then divide P,H,C,lambda by 1-epsilon.',
    physical=list(map(str,P.coeffs())),fourier=list(map(str,H.coeffs())),C_exact=str(C),
    lambda_count_exact=str(lam),U_exact='7802',control_U_exact='7803',threshold_squared_exact='10000/4997',
    normalization_exact=str(H(0)),source_normalization_exact=str(N))
write('candidate_exact.json',variant)
gates={}
for name,p,lo,sign in [('fourier',H,Q(0),1),('physical',P,Q(21,2),-1)]:
    rr=p.numer().complex_roots();assert sum(m for z,m in rr)==p.degree()
    boxes=[]
    for z,m in rr:
        if z.imag.contains(0):
            assert z.imag.is_exact() and z.imag==0
            assert z.real<arb(lo),(name,str(z.real),str(lo))
        boxes.append(dict(real=str(z.real),imag=str(z.imag),multiplicity=m))
    assert sign*p(lo)>0
    gates[name]=dict(pass_gate=True,ray_start_exact=str(lo),required_sign=sign,
        degree=p.degree(),root_multiplicity_sum=sum(m for z,m in rr),all_root_boxes=boxes,
        initial_value_exact=str(p(lo)))
obj=P(0)+C+7802*lam
density=arb(obj)*arb.pi()**8/(2**8*math.factorial(8))
assert density<arb(Q(2492121,10**8))
expressions=[]
for count,upper in [(7802,'0.024920988'),(7803,'0.024921003')]:
    objective=P(0)+C+count*lam
    expression=arb(objective)*arb.pi()**8/(2**8*math.factorial(8))
    assert expression<arb(q(upper))
    expressions.append(dict(U=count,objective_exact=str(objective),
        expression_interval=str(expression),strict_outward_target=upper))
write('expression_controls.json',dict(radial_component_only=True,expressions=expressions))
write('RAYS.json',dict(certified_bound=False,count_conditional=True,gates=gates,
    objective_exact=str(obj),density_expression_interval=str(density),
    comparator_headroom_interval=str(arb(Q(2492121,10**8))-density)))
print('REPAIRED RAYS PASS',str(density),flush=True)

fit=json.loads((D/'inputs/fit.json').read_text());par=list(map(q,fit['parameters']))
angular=json.loads((D/'inputs/angular_dim16.json').read_text());B=Q(angular['B_exact'])
assert C>=0 and lam>=0 and min(par)>=0 and B>0
VK=list(map(F,['1.02','1.05','1.1','1.2','1.4','2']))
LK=list(map(F,['1.2','1.4','1.5','1.6','1.65','1.75']))
threshold=F(5000,4997);endpoint=F(7,4)
cuts=sorted(set([F(1),threshold,endpoint]+[x for x in VK+LK if 1<x<endpoint]))
twopi=2*arb.pi();ratio=arb(C/B);alambda=arb(lam)
AP=arb_poly([arb(x) for x in P.coeffs()])
AP1=arb_poly([arb(x) for x in P.derivative().coeffs()])
AP2=arb_poly([arb(x) for x in P.derivative().derivative().coeffs()])

def profile_parts(mid):
    a=par[0];b=par[7];constant=-arb(par[0])
    for knot,c in zip(VK,par[1:7]):
        if mid>knot:a+=c;constant-=arb(c)*af(knot)
    for knot,c in zip(LK,par[8:]):
        if mid>knot:b+=c;constant-=arb(c)*af(knot).log()
    return arb(a),arb(b),constant

def exp_parts(v,abc):
    a,b,c=abc;phi=a*v+b*v.log()+c;xx=twopi*v
    return xx,(xx-phi).exp(),(xx-2*phi).exp(),a+b/v,-b/v**2

# Negative control: the unrepaired source has a rigorous positive pair residual.
badx=Q(6858771052618467,562949953421312)
badv=arb(badx)/twopi
assert badv>arb(Q(7,4)) and badv<2
aa,bb,cc=profile_parts(F(39,20))
badh=(-(aa*badv+bb*badv.log()+cc)).exp()
bad_residual=arb(P_unrepaired(badx))*(-arb(badx)).exp()-arb(Q(source['C'])/B)*(2*badh-badh*badh)
assert arb(badx)>arb.pi()*arb(Q(10000,4997))
assert bad_residual>0
negative_controls=dict(unrepaired_source_pair=dict(passed=True,x_exact=str(badx),
    positive_residual_interval=str(bad_residual)))

def point_value(v,abc,count):
    vv=af(v);xx,t1,t2,pp,p2=exp_parts(vv,abc)
    ff=AP(xx)-ratio*(2*t1-t2)
    fp=twopi*AP1(xx)-ratio*(2*t1*(twopi-pp)-t2*(twopi-2*pp))
    if count:
        term=alambda*xx.exp();ff-=term;fp-=twopi*term
    return ff,fp

def bound(lo,hi,abc,count):
    mid=(lo+hi)/2;rad=(hi-lo)/2;v=arb(af(mid),af(rad))
    xx,t1,t2,pp,p2=exp_parts(v,abc)
    # Centered polynomial composition avoids dependency blowup from monomial Horner evaluation.
    xm=twopi*af(mid);xr=twopi*af(rad)
    shift=AP2(arb_poly([xm,xr]));poly2=shift[0]
    for k in range(1,len(shift)):poly2+=arb(0,abs(shift[k]))
    f2=twopi**2*poly2-ratio*(2*t1*((twopi-pp)**2-p2)-t2*((twopi-2*pp)**2-2*p2))
    if count:f2-=alambda*twopi**2*xx.exp()
    f0,f1=point_value(mid,abc,count)
    # Taylor's theorem with an interval enclosure of the complete second derivative.
    upper=f0+abs(f1)*af(rad)+abs(f2)*af(rad)**2/2
    return upper

accepted=[];nodes=0;maxdepth=0
for lo,hi in zip(cuts,cuts[1:]):
    abc=profile_parts((lo+hi)/2);count=(lo+hi)/2<threshold
    stack=[(lo,hi,0)]
    while stack:
        a,b,depth=stack.pop();nodes+=1
        if nodes>200000:raise RuntimeError('Bounded pair proof node limit exceeded')
        upper=bound(a,b,abc,count)
        if upper<0:
            accepted.append(dict(lo=str(a),hi=str(b),count_active=count,depth=depth,upper=str(upper)))
            maxdepth=max(maxdepth,depth);continue
        if depth>=36:
            write('FAILURE.json',dict(interval=[str(a),str(b)],depth=depth,upper=str(upper)))
            raise RuntimeError('Pair interval remained unresolved at bounded depth')
        m=(a+b)/2;stack.extend([(m,b,depth+1),(a,m,depth+1)])
    print('PAIR PIECE PASS',str(lo),str(hi),'count',count,'leaves',len(accepted),flush=True)

# Exact rational partition coverage, with half-open count discontinuity handled by both limits.
leaves=sorted(accepted,key=lambda z:F(z['lo']))
def check_partition(cover):
    assert F(cover[0]['lo'])==1 and F(cover[-1]['hi'])==endpoint
    assert all(F(a['hi'])==F(b['lo']) for a,b in zip(cover,cover[1:]))
    assert sum((F(z['hi'])-F(z['lo']) for z in cover),F(0))==endpoint-1
    assert any(F(z['lo'])==threshold and not z['count_active'] for z in cover)
    assert any(F(z['hi'])==threshold and z['count_active'] for z in cover)
    for z in cover:
        a,b=F(z['lo']),F(z['hi'])
        assert a<b and not a<threshold<b
        assert z['count_active']==(b<=threshold)
check_partition(leaves)
for name,mutated in [('deleted_first_interval',leaves[1:]),
    ('incorrect_closed_count_endpoint',[dict(z,count_active=True) if F(z['lo'])==threshold else z for z in leaves])]:
    try:check_partition(mutated)
    except AssertionError:negative_controls[name]=dict(passed=True,corruption_rejected=True)
    else:raise RuntimeError('Negative partition control was not rejected')
write('negative_controls.json',negative_controls)
assert twopi*af(endpoint)>arb(Q(21,2))
# For v>=7/4 all physical values are negative by the stronger exact ray gate.
# h lies in (0,1] for v>=1, so C/B*(2h-h^2) and count allowance are nonnegative.
write('partition.json',leaves)
result=dict(status='PASS',certified_bound=False,count_conditional=True,
    scope='Full conditional radial pair inequality and both rays only. Spherical average-count theorem remains external.',
    finite_pair_domain_v=['1','7/4'],finite_pair_domain_X=['2*pi','7*pi/2'],
    stronger_physical_negative_ray_X='21/2',tail_overlap_proved=True,
    profile_kinks_below_endpoint=list(map(str,sorted(set([x for x in VK+LK if 1<x<=endpoint])))),
    later_kinks_covered_by_physical_negative_ray=['2'],
    count_threshold_v_exact=str(threshold),half_open_count_handling='Left closure proved with count; right closure proved without count, including the threshold.',
    count_U=7802,control_count_U=7803,
    nodes=nodes,accepted_intervals=len(leaves),max_depth=maxdepth,exact_total_width=str(endpoint-1),
    density_expression_interval=str(density),conservative_expression_upper='0.024920988',
    comparator_headroom_interval=str(arb(Q(2492121,10**8))-density),
    
    source_hashes={str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [D/'conditional_refit_exact.json',D/'inputs/fit.json',D/'inputs/angular_dim16.json']})
assert density<arb(q(result['conservative_expression_upper']))
write('RESULT.json',result)
print('FULL CONDITIONAL RADIAL PASS',json.dumps(result),flush=True)
