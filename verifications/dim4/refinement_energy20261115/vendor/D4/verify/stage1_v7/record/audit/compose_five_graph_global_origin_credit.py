#!/usr/bin/env python3
"""Compose exact origin credit with the isolated schema-V6 five-graph proof.

The default mode requires a completed, source-bound global VERIFIED.json,
including every unique closed cyclic branch.  --prepare establishes only
the new local algebra and gives a conditional objective.  Neither mode
replays or replaces the inherited global sign proof.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,importlib.util,itertools,json,math
import certify_origin_credit as radial

ROOT=Path(__file__).resolve().parent
CURRENT=ROOT.parent/'current'
PARENT_COMPOSER=ROOT/'compose_global_origin_credit.py'
PARENT_COMPOSER_SHA256='235e66bd4ca76bf2177134375c1f62ac887a62119592bdd6ab1a190276a198a4'
# The isolated coordinator has a complete independent source and input audit.
COORDINATOR_SHA256='0babb579b3fc84e037c4d3bb45a18f296f633af5afcab535815fcd5ef81d00c0'
TRUSTED_COORDINATORS={COORDINATOR_SHA256:
    CURRENT/'kernel_integration_route/verify_combined_binary_five_antipodal_graphs_small_cap_bounds.py'}
TRUSTED_GRAPH_FAMILIES={COORDINATOR_SHA256:'five'}
TAIL_COORDINATORS=set()
SMALL_CAP_COORDINATORS={COORDINATOR_SHA256}
SMALL_CAP_HELPER=CURRENT/'origin_route/interval_combined_binary_five_graphs_correlated_small_cap_bounds.py'
SMALL_CAP_HELPER_SHA256='3af2b297fb0e13455c42b57cf2f33edba9e22c2dd6c2129f68f51b29aa9890c8'
FIVE_GRAPH_FOLDER=CURRENT/'five_antipodal_graph_route'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def exact(value):
    assert isinstance(value,(str,int)) and not isinstance(value,bool)
    return Q(value)
def dot(a,b):
    assert len(a)==len(b)
    return sum((x*y for x,y in zip(a,b)),Q(0))
def product(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return radial.trim(out)

def verify_radial_nonnegativity(a,kernel,credit):
    degree=kernel['local_degree'];nc=degree+1
    breaks=list(map(exact,kernel['breaks']))
    intervals=[tuple(map(exact,p)) for p in kernel['radial_basis_intervals']]
    assert len(a)==nc*(len(breaks)-1) and len(intervals)==len(breaks)-1
    assert all(x<y for x,y in zip(breaks,breaks[1:]))
    proofs=[]
    for region,(lo,hi) in enumerate(intervals):
        left,right=breaks[region:region+2]
        assert lo<=left<right<=hi
        b=a[nc*region:nc*(region+1)]
        p=radial.bernstein_power(b)
        l,r=(left-lo)/(hi-lo),(right-lo)/(hi-lo)
        detail=dict(region=region,basis_interval=list(map(str,(lo,hi))),
                    physical_closed_interval=list(map(str,(left,right))),
                    normalized_closed_interval=list(map(str,(l,r))),
                    exact_power_coefficients=list(map(str,p)))
        if all(x>=0 for x in b):
            detail['proof']='EXACT_NONNEGATIVE_BERNSTEIN_COEFFICIENTS'
        else:
            sturm=radial.sturm_positive(p,l,r)
            if sturm is not None:
                detail.update(proof='FRESH_EXACT_STURM_STRICT_POSITIVITY',sturm=sturm)
            else:
                # Support the earlier exact factored fallback as well.
                # The supplied factors are merely candidates: their product
                # must equal this actual physical polynomial coefficientwise.
                accepted=False
                for factor in credit.get('radial_nonnegativity',[]):
                    if factor.get('region')!=region or factor.get('formula')!='(1-t)*(t-c)^2*(a+b*t)':continue
                    aa,bb,c=map(exact,(factor['a'],factor['b'],factor['c']))
                    assert aa>=0 and bb>=0 and 0<=c<=1 and 0<=l<=r<=1
                    factored=product(product([Q(1),Q(-1)],[c*c,-2*c,Q(1)]),[aa,bb])
                    assert factored==radial.trim(p)
                    detail.update(proof='FRESH_EXACT_NONNEGATIVE_FACTORIZATION',
                                  formula='(1-t)*(t-c)^2*(a+b*t)',a=str(aa),b=str(bb),c=str(c))
                    accepted=True;break
                assert accepted,('No exact radial nonnegativity proof',region)
        proofs.append(detail)
    return proofs

def verify_local_credit(source,credit):
    required={'origin.json','scalar.json','kernel.json'}
    assert set(credit['source_sha256'])==required
    for name,digest in credit['source_sha256'].items():assert sha(source/name)==digest
    origin=json.loads((source/'origin.json').read_text())
    scalar=json.loads((source/'scalar.json').read_text())
    kernel=json.loads((source/'kernel.json').read_text())
    assert exact(scalar['alpha'])==2 and exact(scalar['G'][0])==1
    assert type(kernel['local_degree']) is int and kernel['local_degree']>=0
    n=(kernel['local_degree']+1)*(len(kernel['breaks'])-1)
    assert kernel['blocks'][0]['dimension']==n+1
    denominator=kernel['common_denominator']
    assert isinstance(denominator,(int,str))
    denominator=int(denominator);assert denominator>0
    M=[[0]*n for i in range(n)]
    for ray in kernel['blocks'][0]['rays']:
        w=ray['weight'];vector=ray['vector']
        assert type(w) is int and w>0 and len(vector)==n+1
        assert all(type(x) is int for x in vector)
        for i,x in enumerate(vector[1:]):
            for j,y in enumerate(vector[1:]):M[i][j]+=w*x*y
    a=list(map(exact,credit['weight_Bernstein_coefficients']))
    b=list(map(exact,credit['range_solution']))
    u=list(map(exact,origin['u']));v=list(map(exact,origin['v']))
    assert len(a)==len(b)==len(u)==len(v)==n
    mat=lambda x:[sum((Q(m,denominator)*q for m,q in zip(row,x)),Q(0)) for row in M]
    assert mat(u)==v and mat(b)==a
    kappa=exact(origin['kappa']);assert dot(u,v)==kappa
    gain=2*dot(u,a)-dot(a,b);assert gain>0
    unew=[x-y for x,y in zip(u,b)];vnew=[x-y for x,y in zip(v,a)]
    assert mat(unew)==vnew
    knew=dot(unew,vnew);assert knew==kappa-gain>=0
    assert list(map(exact,credit['new_psd_origin_u']))==unew
    assert list(map(exact,credit['new_psd_origin_v']))==vnew
    assert exact(credit['new_psd_kappa'])==knew
    assert exact(credit['origin_credit'])==gain and exact(credit['density_credit'])==gain/8
    assert exact(credit['source_objective'])==exact(origin['objective'])
    assert exact(credit['new_formal_objective'])==exact(origin['objective'])-gain/8
    radial_proofs=verify_radial_nonnegativity(a,kernel,credit)
    result=dict(status='FRESH_EXACT_RADIAL_ORIGIN_CREDIT_REPLAY_PASSED',
                matrix_dimension=n,range_equations=n,shifted_range_equations=n,
                density_normalization='alpha=2, G(0)=1, so density objective=(F(0)+A)/8',
                origin_credit=str(gain),density_credit=str(gain/8),
                new_psd_origin_u=list(map(str,unew)),new_psd_origin_v=list(map(str,vnew)),
                new_psd_kappa=str(knew),radial_weight_coefficients=list(map(str,a)),
                range_solution=list(map(str,b)),radial_nonnegativity=radial_proofs,
                local_identity='Enew(C)=Eold(C)-g=(u-b+z_C)^T M (u-b+z_C)+2*sum_x w(|x|)+unchanged other energies',
                nonorigin_kernel_unchanged=True,
                physical_origin_coupling_unchanged='PSD shift -a is exactly restored by the nonnegative pure-linear energy +2*sum w.',
                all_pair_and_cyclic_expressions_unchanged=True,
                finite_range_and_all_cardinality_valid=True)
    return result,origin,kernel,gain

def graph_constant(source):
    path=source/'graph_energies.json'
    assert path.is_file()
    data=json.loads(path.read_text())
    assert data['schema']=='INDEPENDENT_EXTERIOR_GRAPH_ENERGIES_V6'
    assert data['base_input_sha256']=={name:sha(source/name) for name in ('kernel.json','origin.json','scalar.json')}
    weights=list(map(exact,data['weights']))
    assert len(weights)==3 and min(weights)>=0
    near=exact(data['outer_graph_weight']);narrow=exact(data['narrow_graph_weight'])
    assert near>=0 and narrow>=0
    required={'verify_exterior_graphs.py','exterior_graph_certificates.json','graph_repair_bounds.py',
              'verify_outer_exterior_graph.py','outer_exterior_graph_certificate.json','outer_exterior_graph_bounds.py',
              'verify_antipodal_wide_exterior_graph.py','antipodal_wide_exterior_graph_certificate.json',
              'verify_narrow_exterior_graph.py','narrow_exterior_graph_certificate.json','narrow_exterior_graph_bounds.py'}
    assert set(data['source_sha256'])==required
    paths={}
    for name,digest in data['source_sha256'].items():
        assert Path(name).name==name
        path=FIVE_GRAPH_FOLDER/name
        assert sha(path)==digest
        paths[str(path.resolve())]=digest
    constant=dot(weights,list(map(Q,(862,1102,1358))))+120*near+604*narrow
    return constant,data,paths,'five',862

def verify_tail_addition(source,family):
    assert family=='five'
    assert not (source/'gaussian_tail_repair.json').exists(), 'No Gaussian repair is supported by this five-graph coordinator.'
    return None,Q(0),{}

def verify_fifth_graph_freshly(graph):
    """Rebuild the complete E4 count, projection, and exterior proof now."""
    import sys
    common_path=FIVE_GRAPH_FOLDER/'verify_exterior_graphs.py'
    narrow_path=FIVE_GRAPH_FOLDER/'verify_narrow_exterior_graph.py'
    certificate_path=FIVE_GRAPH_FOLDER/'narrow_exterior_graph_certificate.json'
    for path in (common_path,narrow_path,certificate_path):
        assert sha(path)==graph['source_sha256'][path.name]
    previous=sys.modules.get('verify_exterior_graphs')
    previous_write=sys.dont_write_bytecode
    sys.dont_write_bytecode=True
    try:
        spec=importlib.util.spec_from_file_location('verify_exterior_graphs',common_path)
        common=importlib.util.module_from_spec(spec);sys.modules[spec.name]=common
        spec.loader.exec_module(common)
        spec=importlib.util.spec_from_file_location('fresh_composer_fifth_graph',narrow_path)
        narrow=importlib.util.module_from_spec(spec);spec.loader.exec_module(narrow)
        assert narrow.parse is common.parse and narrow.harmonic_power is common.harmonic_power
        assert narrow.prove_nonnegative is common.prove_nonnegative
        proof=narrow.verify(json.loads(certificate_path.read_text()))
    finally:
        sys.dont_write_bytecode=previous_write
        if previous is None:sys.modules.pop('verify_exterior_graphs',None)
        else:sys.modules['verify_exterior_graphs']=previous
    assert proof['status']=='EXACT_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS'
    assert proof['even_ordered_edge_bound']==604
    assert tuple(map(exact,proof['radial_interval']))==(Q(33,28),Q(7,4))
    assert tuple(map(exact,proof['target_angular_interval']))==(Q(-1),Q(-183,200))
    assert proof['count_proof']['count_bound']==151 and proof['degree_proof']['degree_bound']==4
    assert exact(proof['degree_proof']['projected_pair_inner_product_upper'])==-Q(961,1960000)
    assert exact(proof['minimum_selected_distance_squared'])==Q(417087,78400)>Q(7,4)**2
    assert proof['full_triangle_selected_anchor_count']==0
    assert proof['exterior_selected_anchor_count_upper']==1
    assert proof['kernel_diagonal']=='0' and proof['linear_term']=='0'
    assert proof['inherited_neighbor_count_used'] is False
    return proof

def walk_dicts(data):
    if isinstance(data,dict):
        yield data
        for value in data.values():yield from walk_dicts(value)
    elif isinstance(data,list):
        for value in data:yield from walk_dicts(value)

def verify_small_cap_screen(screen,origin):
    assert screen['status']=='EXACT_WHOLE_BRANCH_ZERO_AND_SMALL_CAP_ENCLOSURE_SCREEN'
    caps=origin['pointed_caps'];n=len(caps)
    assert screen['original_count']==n
    zero=screen['zero_cap_indices'];small=screen['small_cap_indices']
    retained=screen['retained_cap_indices']
    combined=zero+small+retained
    assert all(type(i) is int and 0<=i<n for i in combined)
    assert len(combined)==len(set(combined))==n and set(combined)==set(range(n))
    assert screen['retained_count']==len(retained)
    assert screen['omitted_cap_indices']==sorted(zero+small)
    threshold=exact(screen['small_cap_threshold']);assert threshold==Q(1,10**11)
    rows=screen['small_cap_bounds'];assert [r['index'] for r in rows]==small
    total=Q(0)
    for row in rows:
        cap=caps[row['index']]
        magnitude=exact(cap['weight'])*max(map(exact,cap['w']))
        assert exact(row['weight_times_max_radial_coefficient'])==magnitude
        assert 0<=magnitude<threshold
        lower,upper=map(exact,row['whole_branch_cap_bounds'])
        assert 0<=lower<=upper and upper>0
        assert list(map(exact,row['substituted_interval']))==[Q(0),upper]
        total+=upper
    for i in retained:
        cap=caps[i]
        assert exact(cap['weight'])*max(map(exact,cap['w']))>=threshold
    assert exact(screen['small_cap_upper'])==total
    assert list(map(exact,screen['substituted_small_cap_sum_interval']))==[Q(0),total]
    assert screen['method']=='KEEP_ALL_ORIGINAL_VALIDATION; EXACT_ZERO_SCREEN; TINY_NONZERO_CAPS_USE_FIXED_WHOLE_CLOSED_BRANCH_INTERVAL_0_TO_U_IN_EVERY_DESCENDANT_AND_SINGLETON; ADD_SUM_U_TO_ALL_CORRECTION_UPPERS'
    return total

def verify_completed_global(source,origin,kernel,graph,constant,family,capacity0,tail=None,narrow_proof=None):
    assert family=='five' and tail is None and narrow_proof is not None
    path=source/'VERIFIED.json'
    assert path.is_file(),'A completed global VERIFIED.json is required; use --prepare for conditional algebra only.'
    report=json.loads(path.read_text())
    has_tail=tail is not None
    assert report['status']==('EXACT_GLOBAL_BINARY_COMBINED_SIGNED_SPHERICAL_DUAL_GRAPH_GAUSSIAN_TAIL_CERTIFICATE_VERIFIED'
                              if has_tail else 'EXACT_GLOBAL_BINARY_COMBINED_SIGNED_SPHERICAL_GRAPH_COUNT_GAUGE_CERTIFICATE_VERIFIED')
    verifier_digest=report['verifier_sha256']
    assert verifier_digest in TRUSTED_COORDINATORS,'Global coordinator has not been included in this audited transfer contract.'
    assert has_tail==(verifier_digest in TAIL_COORDINATORS)
    verifier=TRUSTED_COORDINATORS[verifier_digest]
    small_caps=verifier_digest in SMALL_CAP_COORDINATORS
    assert sha(verifier)==verifier_digest
    assert TRUSTED_GRAPH_FAMILIES[verifier_digest]==family
    if family in ('four','cap9','dual','near','five'):assert not report['parallel'].get('unproved_jobs')
    expected_inputs={'origin.json','scalar.json','kernel.json','graph_energies.json'}
    if has_tail:expected_inputs.add('gaussian_tail_repair.json')
    expected_codes={'certificate_core.py','verify_scalar.py','audit_identities.py','exact_bernstein_fast.py'}
    assert set(report['input_sha256'])==expected_inputs
    assert set(report['initial_input_sha256'])==expected_inputs|expected_codes
    assert report['initial_input_sha256']==report['final_input_sha256']
    bindings={str(path.resolve()):sha(path)}
    for name,digest in report['initial_input_sha256'].items():
        assert sha(source/name)==digest
        if name in expected_inputs:assert report['input_sha256'][name]==digest
        bindings[str((source/name).resolve())]=digest
    assert report['contract_sha256']==report['final_contract_sha256']
    assert report['contract_sha256'][str(verifier.resolve())]==verifier_digest
    for name,digest in report['contract_sha256'].items():
        assert sha(name)==digest
        bindings[str(Path(name).resolve())]=digest
    if small_caps:
        assert not has_tail
        assert report['contract_sha256'][str(SMALL_CAP_HELPER.resolve())]==SMALL_CAP_HELPER_SHA256
        assert sha(SMALL_CAP_HELPER)==SMALL_CAP_HELPER_SHA256
    if has_tail:
        assert report['contract_sha256'][str(TAIL_HELPER.resolve())]==TAIL_HELPER_SHA256
        assert report['scalar_tail_repair']==tail
        assert exact(report['scalar_tail_density_charge'])==exact(tail['exact_density_charge'])
        assert report['tail_repaired_supported_pair_increment']=='phi-phi_+=min(phi,0)<=0'
        assert report['scalar']['Fourier_half_line']==tail['base_Fourier_proof']
        assert report['scalar']['spatial_tail']==tail['tail_proof']
        assert report['scalar']['proof_method']=='EXACT_BASE_FOURIER_PLUS_NONNEGATIVE_SECOND_GAUSSIAN_AND_REPAIRED_TAIL_COMPARISON'
    # Rebind geometric and signed spherical certificates that the parent
    # replayed separately from its static verifier-source contract.
    count_dir=CURRENT/'base/new_work/spherical_sdp'
    for item in walk_dicts(origin):
        for name_key,hash_key,base in [
            ('certificate','certificate_sha256',count_dir),
            ('spherical_certificate','spherical_certificate_sha256',count_dir),
            ('source_certificate','source_sha256',source/'spherical_energy_sources')]:
            if name_key not in item or hash_key not in item:continue
            name=item[name_key];digest=item[hash_key]
            assert isinstance(name,str) and Path(name).name==name
            target=base/name;assert sha(target)==digest
            bindings[str(target.resolve())]=digest
    nr=len(kernel['breaks'])-1;degree=kernel['angular_degree']
    jobs={(kind,regions) for kind,size in [('full',3),('outside',2)]
          for regions in itertools.combinations_with_replacement(range(nr),size)}
    parallel=report['parallel']
    declared=[(p['kind'],tuple(p['regions'])) for p in parallel['expected_jobs']]
    assert len(declared)==len(set(declared))==len(jobs) and set(declared)==jobs
    assert parallel['number_of_regions']==nr
    assert parallel['expected_full_branches']==math.comb(nr+2,3)
    assert parallel['expected_outside_branches']==math.comb(nr+1,2)
    assert parallel['expected_identity_checks']==6*len(jobs)
    assert parallel['accepted_jobs']==len(jobs)
    assert parallel['coverage_status']=='ALL_EXPECTED_UNIQUE_BRANCHES_AND_IDENTITIES_VERIFIED'
    assert parallel['parent_gates_before_fork']=='FRESH_DEPENDENCIES_SCALAR_ORIGIN_OBJECTIVE_AND_ALL_BASE_PAIRS_PASSED'
    packets=parallel['branch_packets'];assert len(packets)==len(jobs)
    seen=set();cyclic_map={}
    for packet in packets:
        job=packet['kind'],tuple(packet['regions'])
        assert job in jobs and job not in seen;seen.add(job)
        identity=packet['identities'];cyclic=packet['cyclic']
        assert identity['status']=='VERIFIED' and identity['polynomial_identity_checks']==6
        assert identity['angular_degree']==degree
        assert identity['verification_scope']=='EXACT_ARITHMETIC_CROSSCHECKS_ONLY_NOT_A_GLOBAL_CERTIFICATE'
        assert identity['basis_transform_method']=='EXACT_TRANSPOSE_AXIS_CONTRACTION'
        assert identity['interval_converter_implementation_checked'] is False
        assert cyclic['status']=='CERTIFIED' and not cyclic['failures']
        assert (cyclic['kind'],tuple(cyclic['regions']))==job
        if small_caps:verify_small_cap_screen(cyclic['exact_cap_screen'],origin)
        sid=cyclic['signed_spherical_identity']
        assert sid['status']=='EXACT_SIGNED_SPHERICAL_BRANCH_IDENTITIES_VERIFIED'
        assert (sid['kind'],tuple(sid['regions']))==job and sid['checks']==5
        enclosure=cyclic['bernstein_enclosure']
        assert enclosure['method']=='EXACT_OUTWARD_TWO_PRODUCT_BERNSTEIN_V3'
        assert enclosure['all_coefficients_enclosed'] is True and enclosure['coefficient_precision']==192
        exponent=enclosure['output_exponent'];assert type(exponent) is int and -4096<=exponent<=50
        assert exact(enclosure['exact_positive_output_scale'])==Q(2)**exponent>0
        if 'correlated_factor_acceptance' in report['runtime_versions']:
            assert type(cyclic['correlated_factor_leaves']) is int
            assert 0<=cyclic['correlated_factor_leaves']<=cyclic['correction_leaves']<=cyclic['leaves']
        assert list(map(exact,cyclic['graph_weights']))==list(map(exact,graph['weights']))
        assert exact(cyclic.get('outer_graph_weight','0'))==exact(graph.get('outer_graph_weight','0'))
        assert exact(cyclic['narrow_graph_weight'])==exact(graph['narrow_graph_weight'])
        cyclic_map[job]=cyclic
    assert seen==jobs
    assert len(report['cyclic'])==len(jobs)
    assert len({(x['kind'],tuple(x['regions'])) for x in report['cyclic']})==len(jobs)
    assert all(cyclic_map[(x['kind'],tuple(x['regions']))]==x for x in report['cyclic'])
    identity=report['identities'];assert identity['status']=='VERIFIED'
    assert identity['polynomial_identity_checks']==6*len(jobs) and identity['branch_count']==len(jobs)
    assert identity['angular_degree']==degree
    assert report['parent_scalar_and_dependencies']=='FRESH_EXACT_REPLAY_PASSED'
    assert report['origin_identity']=='EXACT_MU_EQUALS_V_AND_VU_EQUALS_KAPPA'
    assert report['full_pair_gauge_nonnegativity']=='ESTABLISHED_BY_THIS_RUNS_FRESH_EXACT_COMPLETE_BASE_PAIR_PROOFS'
    breaks=list(map(exact,kernel['breaks']))
    pairs=report['pairs'];assert pairs and all(p['status']=='CERTIFIED' for p in pairs)
    assert report['complete_pair_majorant_piece_count']==len(pairs)
    # Pair closures must tile every radial branch with no gap or duplicate.
    for region in range(nr):
        intervals=sorted((exact(p['left']),exact(p['right'])) for p in pairs if p['region']==region)
        assert intervals and intervals[0][0]==breaks[region] and intervals[-1][1]==breaks[region+1]
        assert all(a<b for a,b in intervals)
        assert all(b==c for (a,b),(c,d) in zip(intervals,intervals[1:]))
    assert all(type(p['region']) is int and 0<=p['region']<nr for p in pairs)
    energy=report['independent_graph_local_energies']
    assert energy['status']=='FRESH_EXACT_INDEPENDENT_GRAPH_ENERGIES_VERIFIED'
    assert exact(energy['origin_constant'])==constant
    assert list(map(exact,energy['weights']))==list(map(exact,graph['weights']))
    expected_capacities=[capacity0,1102]+([1358] if family in ('cap9','dual','near','five') else [604] if family=='four' else [])
    assert [p['even_ordered_edge_bound'] for p in energy['proofs']]==expected_capacities
    assert all(p['status']=='EXACT_LOCAL_GRAPH_CERTIFICATE_PASS' for p in energy['proofs'][:2])
    if family=='four':
        e=energy['proofs'][2]
        assert e['status']=='EXACT_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS'
        assert tuple(map(exact,e['radial_interval']))==(Q(33,28),Q(7,4))
        assert tuple(map(exact,e['target_angular_interval']))==(Q(-1),Q(-183,200))
        assert e['count_proof']['count_bound']==151 and e['degree_proof']['degree_bound']==4
        assert exact(e['degree_proof']['projected_pair_inner_product_upper'])==-Q(961,1960000)
        assert exact(e['minimum_selected_distance_squared'])==Q(417087,78400)>Q(7,4)**2
    elif family in ('cap9','dual','near','five'):
        e=energy['proofs'][2]
        assert e['status']=='EXACT_PROJECTED_S2_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS'
        assert tuple(map(exact,e['radial_interval']))==(Q(33,28),Q(7,4))
        assert tuple(map(exact,e['target_angular_interval']))==(Q(-1),Q(-7,8))
        assert e['count_proof']['count_bound']==151 and e['degree_proof']['degree_bound']==9
        assert exact(e['degree_proof']['normalized_projected_pair_cosine_upper'])==Q(223,735)
        assert exact(e['degree_proof']['projected_exact_count_upper'])==Q(2560,279)
        assert exact(e['degree_proof']['monotonicity_margin_beta_minus_rho'])==Q(15,392)
        assert exact(e['minimum_selected_distance_squared'])==Q(16335,3136)>Q(7,4)**2
    if 'outer_graph_weight' in graph:
        assert exact(energy['outer_graph_weight'])==exact(graph['outer_graph_weight'])
        outer=energy['outer_graph_proof']
        if family in ('near','five'):
            assert outer['status']=='EXACT_NEAR_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS'
            assert outer['even_ordered_edge_bound']==120
            assert tuple(map(exact,outer['radial_interval']))==(Q(1),Q(103,100))
            assert tuple(map(exact,outer['target_angular_interval']))==(Q(-1),Q(-3,4))
            assert outer['count_proof']['count_bound']==30 and outer['degree_proof']['degree_bound']==4
            assert exact(outer['count_proof']['exact_count_upper'])==Q(288000,9391)
            assert exact(outer['degree_proof']['projected_pair_inner_product_upper'])==-Q(5737,169744)
            assert exact(outer['minimum_selected_distance_squared'])==Q(7,2)>Q(7,4)**2
            assert outer['closed_radial_and_angular_endpoints'] is True
        elif family=='dual':
            assert outer['status']=='EXACT_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS'
            assert outer['even_ordered_edge_bound']==604
            assert tuple(map(exact,outer['radial_interval']))==(Q(33,28),Q(7,4))
            assert tuple(map(exact,outer['target_angular_interval']))==(Q(-1),Q(-183,200))
            assert outer['count_proof']['count_bound']==151 and outer['degree_proof']['degree_bound']==4
            assert exact(outer['degree_proof']['projected_pair_inner_product_upper'])==-Q(961,1960000)
            assert exact(outer['minimum_selected_distance_squared'])==Q(417087,78400)>Q(7,4)**2
        else:
            assert outer['status']=='EXACT_OUTER_EXTERIOR_GRAPH_CERTIFICATE_PASS'
            assert outer['even_ordered_edge_bound']==6252
    assert exact(energy['narrow_graph_weight'])==exact(graph['narrow_graph_weight'])
    assert energy['narrow_graph_proof']==narrow_proof
    assert exact(report['base_objective'])==exact(origin['objective'])
    graph_objective=exact(origin['objective'])+constant/8
    if has_tail:
        assert exact(report['graph_objective'])==graph_objective<exact(graph['threshold'])
        assert exact(report['graph_strict_upper_bound'])==exact(graph['threshold'])
        objective=graph_objective+exact(tail['exact_density_charge'])
        scaled=objective*10**12
        strict_threshold=Q(scaled.numerator//scaled.denominator+1,10**12)
        assert exact(report['strict_upper_bound'])==strict_threshold
    else:
        objective=graph_objective
        assert exact(report['strict_upper_bound'])==exact(graph['threshold'])
    assert exact(report['objective'])==objective<exact(report['strict_upper_bound'])
    return dict(status='COMPLETED_SOURCE_GLOBAL_PROOF_RECORD_AND_BINDING_CHECKED',
                verified_path=str(path),verified_sha256=sha(path),
                verifier_path=str(verifier),verifier_sha256=verifier_digest,
                input_sha256=report['input_sha256'],all_bound_file_sha256=bindings,
                exact_unique_full_branches=math.comb(nr+2,3),
                exact_unique_exterior_branches=math.comb(nr+1,2),
                exact_polynomial_identity_checks=6*len(jobs),
                pair_closed_pieces=len(pairs),source_objective=str(objective),
                source_strict_upper_bound=report['strict_upper_bound'],
                gaussian_tail_proof_freshly_replayed=has_tail,
                sign_proof_replayed=False,
                inheritance='The completed source proof is retained; the freshly proved local identity leaves every pair and cyclic expression exactly unchanged.')

def compose(source,credit_path,output,prepare=False):
    source,credit_path,output=map(lambda x:Path(x).resolve(),(source,credit_path,output))
    assert sha(PARENT_COMPOSER)==PARENT_COMPOSER_SHA256
    code_files=[Path(__file__).resolve(),Path(radial.__file__).resolve(),PARENT_COMPOSER]
    code_hashes={str(p):sha(p) for p in code_files}
    credit_hash=sha(credit_path);credit=json.loads(credit_path.read_text())
    local,origin,kernel,gain=verify_local_credit(source,credit)
    constant,graph,graph_hashes,family,capacity0=graph_constant(source)
    narrow_proof=verify_fifth_graph_freshly(graph)
    tail,tail_charge,tail_bindings=verify_tail_addition(source,family)
    if tail:code_hashes[str(TAIL_HELPER.resolve())]=TAIL_HELPER_SHA256
    source_names=['origin.json','scalar.json','kernel.json','graph_energies.json']
    if tail:source_names.append('gaussian_tail_repair.json')
    source_inputs={str((source/n).resolve()):sha(source/n)
                   for n in source_names}
    inherited=None
    if not prepare or (source/'VERIFIED.json').is_file():
        inherited=verify_completed_global(source,origin,kernel,graph,constant,family,capacity0,tail,narrow_proof)
    old_objective=exact(origin['objective'])+constant/8+tail_charge
    objective=old_objective-gain/8
    threshold=Q((objective*10**12).numerator//(objective*10**12).denominator+1,10**12)
    assert objective<threshold
    proved=inherited is not None and not prepare
    record=dict(status=('EXACT_GLOBAL_ORIGIN_CREDIT_COMPOSITION_VERIFIED' if proved else
                        'EXACT_ORIGIN_CREDIT_COMPOSITION_PREPARED_NOT_GLOBAL'),
                global_density_bound_proved=proved,source=str(source),
                source_input_sha256=source_inputs,
                retained_credit_certificate=str(credit_path),retained_credit_sha256=credit_hash,
                local_credit_replay=local,inherited_global_proof=inherited,
                graph_family=family,graph_capacities=[862,1102,1358,120,604],
                graph_capacity_order=['A','B','E9','near_antipodal','E4'],
                fresh_fifth_graph_proof=narrow_proof,
                fifth_graph_weight=str(exact(graph['narrow_graph_weight'])),
                fifth_graph_density_charge_retained=str(604*exact(graph['narrow_graph_weight'])/8),
                graph_origin_constant_retained=str(constant),graph_density_cost_retained=str(constant/8),
                scalar_tail_repair_replay=tail,scalar_tail_density_charge_retained=str(tail_charge),
                scalar_tail_dependency_sha256=tail_bindings,
                graph_source_sha256=graph_hashes,source_objective=str(old_objective),
                density_credit=str(gain/8),objective=str(objective),objective_decimal=float(objective),
                density_credit_decimal=float(gain/8),
                strict_upper_bound=str(threshold) if proved else None,
                conditional_strict_upper_bound=None if proved else str(threshold),
                threshold_decimal=float(threshold),
                theorem='A nonnegative radial pure-linear local energy restores the PSD-origin shift, reduces the independent origin constant by g, and preserves all physical pair/cyclic expressions.',
                arbitrary_packings=('The unchanged finite-range local-energy transfer theorem applies to arbitrary upper-density sphere packings, inherited from the completed source proof.' if proved else
                                    'Conditional on a completed source global proof: the unchanged finite-range local-energy transfer theorem would apply to arbitrary upper-density sphere packings.'),
                proof_dependency='A conditional preparation is never a global density proof.',
                composition_source_sha256=code_hashes)
    assert sha(credit_path)==credit_hash
    assert code_hashes=={p:sha(p) for p in code_hashes}
    assert source_inputs=={p:sha(p) for p in source_inputs}
    assert graph_hashes=={p:sha(p) for p in graph_hashes}
    assert tail_bindings=={p:sha(p) for p in tail_bindings}
    if inherited:
        assert inherited['all_bound_file_sha256']=={p:sha(p) for p in inherited['all_bound_file_sha256']}
    output.parent.mkdir(parents=True,exist_ok=True)
    temporary=output.with_name(output.name+'.tmp')
    temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(output)
    print(json.dumps({k:record[k] for k in ['status','global_density_bound_proved','objective_decimal',
          'density_credit_decimal','strict_upper_bound','conditional_strict_upper_bound','threshold_decimal']},indent=2))
    return record

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Run without -O.')
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('credit',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--prepare',action='store_true')
    args=p.parse_args();compose(args.source,args.credit,args.output,args.prepare)
