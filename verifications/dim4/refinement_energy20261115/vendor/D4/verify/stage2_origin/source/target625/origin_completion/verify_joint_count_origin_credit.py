#!/usr/bin/env python3
"""Check a stored five-count projection by exact equations and enumeration."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import math
import sys

if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

ROOT=Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upper(n):
    return [(i,j) for i in range(n) for j in range(i,n)]


def verify(data,source):
    assert __debug__
    source=Path(source)
    origin=json.loads((source/'origin.json').read_text())
    kernel=json.loads((source/'kernel.json').read_text())
    assert data['status']=='EXACT_ACCEPTED_KERNEL_COUNT_PROJECTION_AND_HARMONIC_EXTRACTION'
    assert data['global_bound_proved'] is False
    for path,digest in data['source_sha256'].items():assert sha(ROOT/path)==digest
    breaks=list(map(Q,kernel['breaks']))
    assert kernel['dimension']==4 and kernel['local_degree']==6
    assert kernel['radial_basis_intervals']==[['1','3/2']]*6+[['3/2','7/4']]
    assert list(map(Q,data['radii']))==breaks[1:6]
    den=int(kernel['common_denominator']);assert den>0
    M=[[0]*49 for _ in range(49)]
    for ray in kernel['blocks'][0]['rays']:
        assert ray['vector'][0]==0 and type(ray['weight']) is int and ray['weight']>0
        vec=ray['vector'][1:]
        assert len(vec)==49 and all(type(x) is int for x in vec)
        for i in range(49):
            for j in range(49):M[i][j]+=ray['weight']*vec[i]*vec[j]
    A=[[Q(int(region<=cutoff)) for cutoff in range(5)]
       for region in range(7) for coordinate in range(7)]
    B=[list(map(Q,row)) for row in data['indicator_range_solutions']]
    assert len(B)==49 and all(len(row)==5 for row in B)
    for i in range(49):
        for j in range(5):
            assert sum(M[i][k]*B[k][j] for k in range(49))==den*A[i][j]
    u,v=list(map(Q,origin['u'])),list(map(Q,origin['v']))
    assert len(u)==len(v)==49
    assert all(sum(M[i][j]*u[j] for j in range(49))==den*v[i] for i in range(49))
    assert sum(x*y for x,y in zip(u,v))==Q(origin['kappa'])
    D=[[sum(A[k][i]*B[k][j] for k in range(49)) for j in range(5)] for i in range(5)]
    assert D==[list(map(Q,row)) for row in data['D']]
    assert all(D[i][j]==D[j][i] for i in range(5) for j in range(5))
    # Exact LDL positivity; no eigensolver or producer rref result is trusted.
    L=[[Q(int(i==j)) for j in range(5)] for i in range(5)];diagonal=[]
    for i in range(5):
        pivot=D[i][i]-sum(L[i][k]**2*diagonal[k] for k in range(i))
        assert pivot>0
        diagonal.append(pivot)
        for j in range(i+1,5):
            L[j][i]=(D[j][i]-sum(L[j][k]*L[i][k]*diagonal[k] for k in range(i)))/pivot
    G=[list(map(Q,row)) for row in data['D_inverse']]
    assert len(G)==5 and all(len(row)==5 for row in G)
    assert all(sum(D[i][k]*G[k][j] for k in range(5))==int(i==j)
               for i in range(5) for j in range(5))
    m=[sum(A[i][j]*u[i] for i in range(49)) for j in range(5)]
    assert m==list(map(Q,data['count_offsets_m']))
    linear=[sum(G[i][j]*m[j] for j in range(5)) for i in range(5)]
    coefficients=[sum(m[i]*linear[i] for i in range(5))]+linear+[G[i][j] for i,j in upper(5)]
    assert coefficients==list(map(Q,data['projection_quadratic_coefficients']))
    count_records=[origin['nested_count_quadratic']]+origin['conditional_count_polynomials']
    assert len(count_records)==3
    for index,record in enumerate(count_records):
        positions=[list(breaks[1:6]).index(Q(r)) for r in record['radii']]
        assert positions==([2,3,4] if index==0 else [0,1,2,3,4])
        assert record['bounds']==([24,25,26] if index==0 else [24,24,24,25,26])
        if index:
            assert [(x['inner_index'],x['outer_index'],x['maximum']) for x in record['implications']]==[(0,4,24),(1,3,24)]
            assert Q(record['rounding_origin_repair'])==0
        weight=Q(1) if index==0 else Q(record['weight'])
        assert weight>=0
        cs=[weight*Q(x) for x in record['coefficients']]
        n=len(positions);assert len(cs)==1+n+len(upper(n))
        coefficients[0]+=cs[0]
        for i,j in enumerate(positions):coefficients[j+1]+=cs[i+1]
        for coefficient,(i,j) in zip(cs[n+1:],upper(n)):
            coefficients[6+upper(5).index((positions[i],positions[j]))]+=coefficient
    assert coefficients==list(map(Q,data['projection_plus_existing_count_coefficients']))
    denominator=math.lcm(*(x.denominator for x in coefficients))
    ints=[int(x*denominator) for x in coefficients]
    # Independent explicit nested loops, evaluating the augmented quadratic.
    H=[[0]*6 for _ in range(6)];H[0][0]=ints[0]
    for i in range(5):H[0][i+1]=H[i+1][0]=ints[i+1]
    for value,(i,j) in zip(ints[6:],upper(5)):H[i+1][j+1]=H[j+1][i+1]=value
    profiles={};checked=0
    for n4 in range(27):
        for n3 in range(min(n4,25)+1):
            for n2 in range(min(n3,24)+1):
                for n1 in range(n2+1):
                    if n1==24 and n3!=24:continue
                    for n0 in range(n1+1):
                        if n0==24 and n4!=24:continue
                        counts=(n0,n1,n2,n3,n4);x=(1,)+counts
                        value=sum(x[i]*(H[i][i]*x[i]+2*sum(H[i][j]*x[j] for j in range(i+1,6))) for i in range(6))
                        checked+=1
                        if n4 not in profiles or value<profiles[n4][0]:profiles[n4]=(value,counts)
    assert checked==data['count_vectors_checked']==165503
    stored=data['conditional_outer_count_profile'];assert len(stored)==27
    assert [x['N'] for x in stored]==list(range(27))
    for row in stored:
        value,counts=profiles[row['N']]
        assert Q(value,denominator)==Q(row['exact_floor'])
        # A stored minimizer must attain the same value; tied minima are harmless.
        target=(1,)+tuple(row['minimizing_counts'])
        assert tuple(sorted(target[1:]))==target[1:] and target[-1]==row['N']
        assert sum(target[i]*sum(H[i][j]*target[j] for j in range(6)) for i in range(6))==value
    floor=min(value for value,counts in profiles.values())
    credit=Q(floor,denominator)
    assert credit==Q(data['universal_projection_plus_count_floor']) and credit>0
    return dict(status='EXACT_FIVE_COUNT_ORIGIN_CREDIT_LOCAL_IDENTITY_AND_FLOOR_VERIFIED',
                exact_origin_credit=str(credit),exact_density_credit=str(credit/8),
                radial_indicator_identity='Each of the first five cumulative seven-coordinate Bernstein blocks sums to its exact closed radial count indicator.',
                radial_endpoint_convention='First piece [1,b1]; later pieces (bi,b_(i+1)]. The inherited global proof also covers every opposite one-sided closure.',
                residual_square_identity='E0-(N+m)^T*D_inverse*(N+m)=(z+u-B*D_inverse*(N+m))^T*M*(z+u-B*D_inverse*(N+m))>=0',
                indicator_range_equations=245,origin_range_equations=49,inverse_equations=25,
                exact_positive_LDL_pivots=list(map(str,diagonal)),count_vectors_checked=checked,
                conditional_outer_count_profile=[dict(N=N,exact_floor=str(Q(value,denominator)),minimizer=list(counts))
                                                 for N,(value,counts) in sorted(profiles.items())],
                exact_total_nonnegative_energy='Original E0 plus the three unchanged count energies minus the new origin credit.',
                all_physical_origin_couplings_unchanged=True,all_nonorigin_kernel_entries_unchanged=True,
                all_pair_and_cyclic_expressions_unchanged=True,
                cardinality_geometry_must_be_freshly_verified=True,global_bound_proved=False,
                verifier_sha256=sha(Path(__file__)))


if __name__=='__main__':
    source=ROOT/'current/candidates/round9_six_graphs_degree18'
    path=Path(__file__).with_name('ACCEPTED_JOINT_COUNT_PROJECTION.json')
    result=verify(json.loads(path.read_text()),source)
    result['projection_certificate_sha256']=sha(path)
    output=path.with_name('LOCAL_JOINT_COUNT_CREDIT_VERIFIED.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],credit=float(Q(result['exact_origin_credit'])),
                         output_sha256=sha(output)),indent=2))
