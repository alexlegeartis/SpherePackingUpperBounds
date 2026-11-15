"""Exact local graph certificates; never a global sphere-packing acceptance.

Run with Python 3.12+ (no third-party dependencies). Every off-diagonal
angular interval is proved by rational Bernstein subdivision. Harmonic
matrices are reconstructed from explicit positive rank-one rays.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import hashlib,json,sys

def add(a,b):
    c=[Q(0)]*max(len(a),len(b))
    for k,v in enumerate(a):c[k]+=v
    for k,v in enumerate(b):c[k]+=v
    return c

def zonals(degree):
    u=[[Q(1)],[Q(0),Q(2)]]
    for k in range(2,degree+1):u.append(add([Q(0)]+[2*v for v in u[-1]],[-v for v in u[-2]]))
    return [[v/Q(k+1) for v in z] for k,z in enumerate(u[:degree+1])]

def harmonic_power(a):
    p=[Q(0)]
    for c,z in zip(a,zonals(len(a)-1)):p=add(p,[c*v for v in z])
    return p

def bernstein(p,left,right):
    n=len(p)-1;t=[Q(0)]*(n+1)
    for j,c in enumerate(p):
        for k in range(j+1):t[k]+=c*comb(j,k)*left**(j-k)*(right-left)**k
    return [sum((t[k]*Q(comb(i,k),comb(n,k)) for k in range(i+1)),Q(0)) for i in range(n+1)]

def halves(b):
    rows=[b]
    while len(rows[-1])>1:rows.append([(u+v)/2 for u,v in zip(rows[-1],rows[-1][1:])])
    return [r[0] for r in rows],[r[-1] for r in reversed(rows)]

def prove_nonnegative(p,left,right):
    assert left<=right
    stack=[(bernstein(p,left,right),0)];nodes=leaves=depth_max=0
    while stack:
        b,d=stack.pop();nodes+=1;depth_max=max(d,depth_max)
        assert nodes<=200000 and d<80,'incomplete sign proof'
        if min(b)>=0:leaves+=1;continue
        assert b[0]>=0 and b[-1]>=0,'negative exact endpoint'
        assert max(b)>=0,'negative whole interval'
        l,r=halves(b);stack.extend([(r,d+1),(l,d+1)])
    return dict(nodes=nodes,leaves=leaves,depth=depth_max)

def parse(s):
    assert isinstance(s,str)
    q=Q(s);assert str(q)==s
    return q

def max_angle(a,b):
    # For fixed s, the derivative numerator r^2-s^2+1 increases in r;
    # any interior critical point is a minimum. Iterate in both variables.
    return max((r*r+s*s-1)/(2*r*s) for r in a for s in b)

def verify(d):
    assert d['schema']=='EXTERIOR_TWO_COLOR_GRAPH_V1'
    assert d['global_density_bound_proved'] is False
    intervals=[list(map(parse,s)) for s in d['radial_closure_intervals']]
    a,b=intervals
    assert Q(1)<=a[0]<=a[1]<=b[0]<=b[1]<=Q(7,4)
    left,right=map(parse,d['target_angular_interval'])
    assert -1<=left<right<0
    # On this negative angular band, distance^2 is increasing in both radii.
    min_distance_squared=a[0]**2+b[0]**2-2*a[0]*b[0]*right
    assert min_distance_squared>Q(49,16)
    cutoffs=[max_angle(a,a),max_angle(a,b),max_angle(b,b)]
    assert right<=cutoffs[1] and max(cutoffs)<1
    vectors=[list(map(parse,v)) for v in d['psd_ray_vectors']]
    assert all(len(v)==2 for v in vectors)
    weights=[list(map(parse,w)) for w in d['psd_ray_weights']]
    assert weights and all(len(w)==len(vectors) and min(w)>=0 for w in weights)
    matrices=[list(map(parse,d['constant_matrix_upper_triangle']))]
    for w in weights:
        matrices.append([sum((c*v[i]*v[j] for c,v in zip(w,vectors)),Q(0)) for i,j in ((0,0),(0,1),(1,1))])
    assert matrices==[list(map(parse,m)) for m in d['harmonic_matrices']]
    polys=[harmonic_power([m[k] for m in matrices]) for k in range(3)]
    domains=[(0,Q(-1),cutoffs[0],0),(2,Q(-1),cutoffs[2],0),(1,Q(-1),left,0),(1,left,right,1),(1,right,cutoffs[1],0)]
    signs=[]
    for k,lo,hi,target in domains:
        p=[-v for v in polys[k]];p[0]-=target
        signs.append(dict(entry=k,left=str(lo),right=str(hi),target=target,proof=prove_nonnegative(p,lo,hi)))
    aa,bb,cc=matrices[0];det=aa*cc-bb*bb
    assert aa>0 and cc>0 and det>0
    fN=sum(m[0] for m in matrices);fM=sum(m[2] for m in matrices)
    ncenter=(cc*fN-bb*fM)/(2*det)
    realcap=(cc*fN*fN-2*bb*fN*fM+aa*fM*fM)/(4*det)
    schur=det/cc
    limit=d['unrestricted_integer_n_cutoff']
    assert isinstance(limit,int) and limit>0
    # For all real m, q(n,m)<=realcap-schur*(n-ncenter)^2.
    # This bound is decreasing for n>=limit and already negative there.
    assert ncenter<limit and realcap-schur*(limit-ncenter)**2<0
    candidates=[]
    for n in range(limit):
        mcenter=(fM-2*bb*n)/(2*cc)
        m0=max(0,mcenter.numerator//mcenter.denominator)
        for m in sorted(set([0,m0,m0+1])):
            value=n*fN+m*fM-aa*n*n-2*bb*n*m-cc*m*m
            candidates.append((value,n,m))
    cap,n,m=max(candidates)
    bound=2*(cap.numerator//(2*cap.denominator))
    assert bound==d['claimed_even_ordered_edge_bound']
    return dict(name=d['name'],status='EXACT_LOCAL_GRAPH_CERTIFICATE_PASS',radial_closure_intervals=d['radial_closure_intervals'],target_angular_interval=d['target_angular_interval'],spherical_cutoffs=list(map(str,cutoffs)),minimum_selected_distance_squared=str(min_distance_squared),whole_interval_sign_proofs=signs,power_coefficients=[[str(c) for c in p] for p in polys],constant_matrix=list(map(str,matrices[0])),constant_matrix_determinant=str(det),diagonal_at_one=[str(fN),str(fM)],real_quadratic_maximum=str(realcap),unrestricted_integer_n_cutoff=limit,integer_candidates_tested=len(candidates),maximizer=[n,m],exact_ordered_edge_upper=str(cap),even_ordered_edge_bound=bound,gap_to_next_even=str(bound+2-cap),inherited_neighbor_count_used=False,candidate_applied=False,global_density_bound_proved=False)

if __name__=='__main__':
    p=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('exterior_graph_certificates.json')
    data=json.loads(p.read_text());reports=[verify(d) for d in data['certificates']]
    print(json.dumps(dict(status='EXACT_EXTERIOR_GRAPH_CERTIFICATES_PASS',global_density_bound_proved=False,certificate_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),certificates=reports),indent=2))
