"""Independent exact local graph proof and isolated V7 helper bookkeeping."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import ast
import copy
import hashlib
import importlib.util
import itertools
import json
import sys

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent.parent
HERE=ROOT/'geometry'
ROUTE=ROOT/'current/six_exterior_graph_route'
CONSTRUCTION=HERE/'SIX_GRAPHS_SMALL_CAP_BOUNDS_HELPER_CONSTRUCTION.json'
OUTPUT=HERE/'six_graphs_small_cap_bounds_helper_exact_audit.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def nested(path):
    tree=ast.parse(path.read_text())
    check=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='check')
    inner=next(n for n in check.body if isinstance(n,ast.FunctionDef) and n.name=='all_corrections')
    return check,inner


def isolated(inner,path,kind,weights,near_weight,narrow_weight,moderate_weight,modules):
    graph,near,narrow,moderate=modules
    namespace=dict(Q=Q,lru_cache=lru_cache,kind=kind,graph_weights=weights,
        outer_graph_weight=near_weight,narrow_graph_weight=narrow_weight,moderate_graph_weight=moderate_weight,
        cyclic_decrement_lower=graph.cyclic_decrement_lower,
        outer_graph_lower=near.cyclic_decrement_lower,
        narrow_graph_lower=narrow.cyclic_decrement_lower,
        moderate_graph_lower=moderate.cyclic_decrement_lower,
        regions=(5,5),tents=(),caps=(),core=None,basis_intervals=None,probes=(),budget=None,
        allocations=(),small_cap_upper=Q(1,10**17),
        tent_bounds=lambda *args:(Q(2,17),Q(3,17)),
        cap_bounds=lambda *args:(Q(0),Q(0)),
        lift_bounds=lambda *args:(Q(5,19),Q(7,19)),
        allocation_bounds=lambda *args:(Q(11,23),Q(13,23)))
    tree=ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(inner)],type_ignores=[]))
    exec(compile(tree,str(path),'exec'),namespace)
    return namespace['all_corrections']


def coordinate(r,s,t):return (r,s,(1+t)*2*r*s/((r+s)**2-Q(49,16)))


def indicator(point):
    r,s,z=point;t=-1+z*((r+s)**2-Q(49,16))/(2*r*s)
    return int(Q(7,5)<=r<=Q(3,2) and Q(7,5)<=s<=Q(3,2) and Q(-4,5)<=t<=Q(-11,20))


def audit():
    construction=json.loads(CONSTRUCTION.read_text())
    source=Path(construction['original_path']);target=Path(construction['new_path'])
    assert sha(source)==construction['original_sha256']=='3af2b297fb0e13455c42b57cf2f33edba9e22c2dd6c2129f68f51b29aa9890c8'
    assert sha(target)==construction['new_sha256']=='0c34525decc059471e69a3d33361a0f57bbe9503a4cdffc7d670b3ad89b707f5'
    reverse=target.read_text()
    for edit in reversed(construction['changes']):
        assert reverse.count(edit['new'])==1
        reverse=reverse.replace(edit['new'],edit['old'])
    assert reverse==source.read_text()
    scalar_path=ROOT/'current/candidates/round9_six_graphs_degree18/verify_scalar.py'
    paths=[source,target,CONSTRUCTION,scalar_path,*sorted(ROUTE.glob('*.py')),*sorted(ROUTE.glob('*.json'))]
    bindings={str(path):sha(path) for path in paths}
    sys.path.insert(0,str(ROUTE))
    graph=load('graph_repair_bounds',ROUTE/'graph_repair_bounds.py')
    near=load('outer_exterior_graph_bounds',ROUTE/'outer_exterior_graph_bounds.py')
    narrow=load('narrow_exterior_graph_bounds',ROUTE/'narrow_exterior_graph_bounds.py')
    moderate=load('moderate_exterior_graph_bounds',ROUTE/'moderate_exterior_graph_bounds.py')
    verifier=load('_independent_moderate_local_verifier',ROUTE/'verify_moderate_exterior_graph.py')
    proof=verifier.verify(json.loads((ROUTE/'moderate_exterior_graph_certificate.json').read_text()))
    assert proof['status']=='EXACT_SINGLE_COLOR_EXTERIOR_GRAPH_CERTIFICATE_PASS'
    assert proof['even_ordered_edge_bound']==moderate.CAPACITY==1880
    assert tuple(map(Q,proof['radial_interval']))==moderate.RADIAL==(Q(7,5),Q(3,2))
    assert tuple(map(Q,proof['target_angular_interval']))==moderate.ANGULAR==(Q(-4,5),Q(-11,20))
    assert Q(proof['minimum_selected_distance_squared'])==Q(1519,250)>moderate.KERNEL_SUPPORT_RADIUS**2
    # Independently prove the same strict interval signs by primitive
    # pseudo-remainder Sturm sequences, without the Bernstein checker.
    sturm=load('_independent_scalar_sturm_for_graph',scalar_path)
    power=list(map(Q,proof['power_coefficients']));sturm_proofs=[]
    for row in proof['whole_interval_sign_proofs']:
        left,right=Q(row['left']),Q(row['right'])
        polynomial=[-x for x in power];polynomial[0]-=row['target']
        sequence=sturm.sturm(sturm.integers(polynomial))
        lsigns=[sturm.sign_at(p,left) for p in sequence]
        rsigns=[sturm.sign_at(p,right) for p in sequence]
        count=sturm.variations(lsigns)-sturm.variations(rsigns)
        assert count==0 and lsigns[0]>0 and rsigns[0]>0
        sturm_proofs.append(dict(left=str(left),right=str(right),target=row['target'],
            roots=count,left_signs=lsigns,right_signs=rsigns))
    a0=Q(proof['edge_count_proof']['constant_coefficient'])
    diagonal=Q(proof['edge_count_proof']['diagonal_at_one'])
    q=lambda n:n*diagonal-a0*n*n
    assert diagonal/(2*a0)<400 and q(400)<0
    maximum=max(q(n) for n in range(401))
    assert maximum==Q(proof['edge_count_proof']['exact_integer_cardinality_edge_upper'])
    assert 2*(maximum.numerator//(2*maximum.denominator))==1880
    old_check,old_inner=nested(source);check,inner=nested(target)
    args=[arg.arg for arg in check.args.args]
    defaults=dict(zip(args[-len(check.args.defaults):],check.args.defaults))
    assert ast.unparse(defaults['moderate_graph_weight'])=='Q(0)'
    validation=[node for node in check.body if isinstance(node,ast.Assert) and 'moderate_graph_weight' in ast.unparse(node)]
    assert len(validation)==1
    validation=compile(ast.fix_missing_locations(ast.Module(body=validation,type_ignores=[])),str(target),'exec')
    invalid_controls=[]
    for value in [Q(-1,10**8),0,0.0,False,'0']:
        try:exec(validation,dict(Q=Q,moderate_graph_weight=value))
        except AssertionError:invalid_controls.append(repr(value))
        else:raise AssertionError('Invalid weight accepted')
    for value in [Q(0),Q(1,10**8)]:exec(validation,dict(Q=Q,moderate_graph_weight=value))
    modules=(graph,near,narrow,moderate);eps=Q(1,10**9)
    radii=[Q(1),Q(103,100),Q(7,5)-eps,Q(7,5),Q(1152597,800000),Q(3,2),Q(3,2)+eps,Q(7,4)]
    cosines=[Q(-1),Q(-183,200),Q(-7,8),Q(-4,5)-eps,Q(-4,5),Q(-4,5)+eps,
             Q(-27,40),Q(-11,20)-eps,Q(-11,20),Q(-11,20)+eps,Q(-1,2)]
    singleton_count=selected_count=0;boxes=[]
    for weights,near_weight,narrow_weight,moderate_weight in [
        ((Q(0),Q(0),Q(1,50000000)),Q(1,5000000),Q(1,5000000),Q(1,10**8)),
        ((Q(1,10**7),Q(2,10**7),Q(3,10**7)),Q(4,10**7),Q(5,10**7),Q(6,10**7))]:
        old=isolated(old_inner,source,'outside',weights,near_weight,narrow_weight,Q(0),modules)
        zero=isolated(inner,target,'outside',weights,near_weight,narrow_weight,Q(0),modules)
        new=isolated(inner,target,'outside',weights,near_weight,narrow_weight,moderate_weight,modules)
        for r,s,t in itertools.product(radii,radii,cosines):
            point=coordinate(r,s,t)
            if not 0<=point[2]<=1:continue
            box=tuple((x,x) for x in point)
            expected=moderate_weight*indicator(point)
            assert moderate.cyclic_decrement_lower('outside',box,moderate_weight)==expected
            assert moderate.cyclic_decrement_upper('outside',box,moderate_weight)==expected
            assert zero(box)==old(box)
            assert tuple(a-b for a,b in zip(new(box),old(box)))==(expected,)*4
            singleton_count+=1;selected_count+=int(expected>0)
        r=s=Q(36,25)
        samples=[
            ('strict_selected',((r,r+eps),(s,s+eps),(coordinate(r,s,Q(-27,40))[2],coordinate(r,s,Q(-27,40))[2]+eps))),
            ('lower_angular_straddle',((r,r),(s,s),(coordinate(r,s,Q(-4,5)-eps)[2],coordinate(r,s,Q(-4,5)+eps)[2]))),
            ('upper_angular_straddle',((r,r),(s,s),(coordinate(r,s,Q(-11,20)-eps)[2],coordinate(r,s,Q(-11,20)+eps)[2]))),
            ('lower_radial_straddle',((Q(7,5)-eps,Q(7,5)+eps),(s,s),(coordinate(Q(7,5),s,Q(-27,40))[2],)*2)),
            ('upper_radial_straddle',((Q(3,2)-eps,Q(3,2)+eps),(s,s),(coordinate(Q(3,2),s,Q(-27,40))[2],)*2)),
        ]
        for name,box in samples:
            lo=moderate.cyclic_decrement_lower('outside',box,moderate_weight)
            hi=moderate.cyclic_decrement_upper('outside',box,moderate_weight)
            assert hi==moderate_weight and zero(box)==old(box)
            assert tuple(a-b for a,b in zip(new(box),old(box)))==(lo,hi,lo,lo)
            count=0
            for point in itertools.product(*[(left,(left+right)/2,right) for left,right in box]):
                assert lo<=moderate_weight*indicator(point)<=hi
                count+=1
            boxes.append(dict(name=name,bounds=[[str(a),str(b)] for a,b in box],
                exact_increment_interval=list(map(str,(lo,hi))),exact_direct_point_controls=count))
        old_full=isolated(old_inner,source,'full',weights,near_weight,narrow_weight,Q(0),modules)
        new_full=isolated(inner,target,'full',weights,near_weight,narrow_weight,moderate_weight,modules)
        for box in [((Q(1),Q(7,4)),)*3,((Q(3,2),Q(3,2)),(Q(7,4),Q(7,4)),(Q(7,4),Q(7,4)))]:
            assert old_full(box)==new_full(box)
    failure=(Q(1152597,800000),Q(1152597,800000),Q(33,128))
    failure_decrement=moderate.cyclic_decrement_lower('outside',tuple((x,x) for x in failure),Q(1,10**8))
    assert failure_decrement==Q(1,10**8)
    old_charge=Q(4299,200000000);new_charge=moderate.density_cost(Q(1,10**8))
    assert new_charge==Q(47,20000000)
    assert all(sha(Path(path))==digest for path,digest in bindings.items())
    return dict(status='EXACT_SIXTH_GRAPH_PROOF_AND_HELPER_BOOKKEEPING_AUDIT_PASS',
        source_contract=bindings,exact_reverse_source_reconstruction=True,
        fresh_local_graph_proof=proof,independent_whole_interval_sturm_proofs=sturm_proofs,
        independent_integer_quadratic_count=401,quadratic_negative_decreasing_from_cardinality=400,
        default_zero_exactly_recovers_parent=True,invalid_weight_controls=invalid_controls,
        singleton_controls=singleton_count,selected_singletons=selected_count,box_controls=boxes,
        full_branches_exactly_unchanged=True,existing_cap_enclosure_and_singleton_rules_unchanged=True,
        known_V6_failure_point=list(map(str,failure)),known_V6_failure_graph_decrement=str(failure_decrement),
        new_graph_density_charge=str(new_charge),combined_graph_density_charge=str(old_charge+new_charge),
        exact_point_feasibility_replay_performed=False,
        scope='Independent exact local polynomial/Sturm and integer-cardinality proof; actual all_corrections AST with real graph functions and fixed unrelated intervals isolates additive bookkeeping. No complete sign, scalar, pair or cyclic replay.',
        frozen_sources_or_candidate_inputs_modified=False,global_density_bound_proved=False)


if __name__=='__main__':
    assert not OUTPUT.exists()
    report=audit();OUTPUT.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ('status','singleton_controls','selected_singletons',
        'new_graph_density_charge','combined_graph_density_charge')},indent=2))
