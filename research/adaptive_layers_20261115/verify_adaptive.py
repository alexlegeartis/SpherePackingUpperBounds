"""Exact adaptive counting: at most one exceptional binary coset, then fit M."""
from itertools import product,combinations
from fractions import Fraction as F
from flint import fmpq_mat,fmpz_mat
import json,argparse
from common import P,theta,root_gram,metric_audit,finish
if not __debug__:raise RuntimeError('Run without -O.')


def bitvector(b,r=5):return [(b>>i)&1 for i in range(r)]

def quadratic(b,Q):
    a=bitvector(b,len(Q))
    return sum(a[i]*Q[i][j]*a[j]for i in range(len(Q))for j in range(len(Q)))//2%2


def embedded_model(Q):
    r=len(Q);counts=[0,0]
    U=[[Q[i][i]//2 if i==j else Q[i][j]%2 if i<j else 0 for j in range(r)]for i in range(r)]
    images=set()
    for b in range(1<<r):
        a=bitvector(b,r);y=[sum(U[i][j]*a[j]for j in range(r))%2 for i in range(r)]
        assert sum(x*z for x,z in zip(a,y))%2==quadratic(b,Q)
        images.add((tuple(a),tuple(y)))
        if b:counts[quadratic(b,Q)]+=1
    assert len(images)==1<<r
    return counts


def adaptive77():
    t,_=theta(72);pop=2**71-2**35-1
    bad5=sum(t[j]for j in(4,6,8))//2;bad4=sum(t[j]for j in(4,6))//2
    cost=F(31*(bad5+bad4),pop);assert cost<2
    # Z<2 guarantees all classes have horizontal norm >=4, with at most
    # one below5. A totally singular five-space exists in either ambient type.
    Q=root_gram(5);M=[[F(3,2)*x for x in v]for v in Q];exception=24
    # The class (0,0,0,1,1) is represented physically by 2*e1, not a root.
    targets=[F(8)if b==0 else F(4)if b==exception else F(3)for b in range(32)]
    metric=metric_audit(M,targets)
    assert F(metric['minimum_by_parity'][str(exception)])==6
    # Explicit GF(2) involutions send the exceptional vertical class to any
    # nonzero class of the totally singular horizontal model.
    for bad in range(1,32):
        if bad==exception:transform=lambda x:x
        else:
            mask=next(x for x in range(1,32)if(x&bad).bit_count()%2 and(x&exception).bit_count()%2)
            w=bad^exception
            transform=lambda x:x^(w if(x&mask).bit_count()%2 else 0)
        assert transform(exception)==bad
        assert len({transform(x)for x in range(32)})==32
        assert all(transform(x^y)==transform(x)^transform(y)for x in range(32)for y in range(32))
    d2=F(2**77)/F(metric['determinant']);assert d2==F(2**80,3**5)
    return finish(77,'Gamma72 adaptive singular cosets above D5',d2,'new_lattice_existence_proof',
                  {'parent_dimension':72,'parent_minimum':8,'model_type_counts':[31,0],
                   'weighted_expected_cost':str(cost),'integer_cost_strictly_below':2,
                   'all_classes_horizontal_squared_norm_at_least':'4',
                   'all_but_at_most_one_at_least':'5','exception_parity':exception,
                   'all_31_exception_reparametrizations_checked':True,'metric_audit':metric})


def adaptive85():
    t,_=theta(80);pop=2**79-2**39-1;Q=root_gram(5)
    assert embedded_model(Q)==[11,20]
    bad0=sum(t[j]for j in(4,6,8))//2
    bad1_hi=sum(t[j]for j in(5,7,9))//2
    bad1_lo=sum(t[j]for j in(5,7))//2
    cost=F(22*bad0+20*(bad1_hi+bad1_lo),pop);assert cost<2
    # Cost2 for a singular class below5, cost2 for a nonsingular class
    # below9/2, cost1 for one below11/2. Thus at most one shallow nonsingular.
    M=json.loads((P/'adaptive85_metric.json').read_text())['gram']
    targets=[F(8)if b==0 else F(7,2)if b==1 else F(5,2)if quadratic(b,Q)else F(3)for b in range(32)]
    metric=metric_audit(M,targets);assert F(metric['determinant'])==F(3875,256)
    R=fmpq_mat([[int(j==i)-int(j==i+1)for j in range(5)]for i in range(4)]+[[int(j in(3,4))for j in range(5)]])
    assert(R*R.transpose()).tolist()==fmpq_mat(Q).tolist()
    Ri=R.inv();classes={}
    for i,j in combinations(range(5),2):
        for a,b in product((-1,1),repeat=2):
            v=[0]*5;v[i]=a;v[j]=b;c=fmpq_mat([v])*Ri
            assert all(x.denominator==1 for x in c.tolist()[0])
            bits=sum((int(c[0,k])%2)<<k for k in range(5));classes[bits]=(v,c)
    assert set(classes)=={b for b in range(1,32)if quadratic(b,Q)}and len(classes)==20
    for bits,(v,c)in classes.items():
        i,j=[k for k,x in enumerate(v)if x];others=[k for k in range(5)if k not in(i,j)]
        U=fmpq_mat(5,5);U[0,i]=v[i];U[1,j]=-v[j]
        for row,col in enumerate(others,2):U[row,col]=1
        C=R*U*Ri
        assert all(x.denominator==1 for vv in C.tolist()for x in vv)
        assert abs(C.det())==1 and C*fmpq_mat(Q)*C.transpose()==fmpq_mat(Q)
        assert C.tolist()[0]==c.tolist()[0]
        for old in range(32):
            z=fmpq_mat([bitvector(old)])*C
            new=sum((int(z[0,k])%2)<<k for k in range(5))
            assert quadratic(new,Q)==quadratic(old,Q)
    d2=F(2**85)/F(metric['determinant']);assert d2==F(2**93,3875)
    return finish(85,'L80 adaptive exceptional coset and fitted D5 metric',d2,'new_lattice_existence_proof',
                  {'parent_dimension':80,'parent_minimum':8,'model_type_counts':[11,20],
                   'weighted_expected_cost':str(cost),'integer_cost_strictly_below':2,
                   'singular_horizontal_squared_norm_at_least':'5',
                   'nonsingular_horizontal_squared_norm_at_least':'9/2',
                   'all_but_at_most_one_nonsingular_at_least':'11/2',
                   'all_20_exception_root_classes_and_isometries_checked':True,'metric_audit':metric})


def audit():return {'method':'Weighted integer bad-coset cost below2, followed by an adapted vertical metric',
                    'rows':[adaptive77(),adaptive85()],'explicit_representatives':False,
                    'all_arithmetic_and_infinite_layer_inequalities_exact':True,
                    'inherited_inputs':'Published Gamma72/L80 minimum8, unimodularity and extremal theta series; standard binary quadratic-space transitivity'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args();r=audit()
    if args.write_result:(P/'adaptive_results.json').write_text(json.dumps(r,indent=2)+'\n')
    else:assert json.loads((P/'adaptive_results.json').read_text())==r
    for a in r['rows']:print(a['dimension'],a['center_density'],'gain',a['gain_over_checkpoint_percent'])
