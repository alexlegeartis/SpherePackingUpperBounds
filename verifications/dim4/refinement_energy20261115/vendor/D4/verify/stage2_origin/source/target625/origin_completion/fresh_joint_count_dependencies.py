#!/usr/bin/env python3
"""Freshly replay the accepted three count energies and their geometry."""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'current/base/new_work'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify(source,accepted,output_dir):
    assert __debug__
    source,output_dir=Path(source),Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)
    old=json.loads(Path(accepted).read_text())
    assert old['status']=='EXACT_GLOBAL_MONOTONE_GRAPH_EXTENSION_AND_ORIGIN_CREDIT_VERIFIED'
    assert old['global_density_bound_proved'] is True
    protected=old['proof_file_sha256']
    by_resolved={str((ROOT/p).resolve()):digest for p,digest in protected.items()}
    binding={}

    def bind(path):
        digest=sha(path)
        assert by_resolved[str(path.resolve())]==digest,str(path)
        binding[str(path.relative_to(ROOT))]=digest
        return path

    origin_path=bind(source/'origin.json')
    kernel_path=bind(source/'kernel.json')
    origin=json.loads(origin_path.read_text())
    breaks=list(map(Q,json.loads(kernel_path.read_text())['breaks']))
    spherical=load('joint_count_fresh_spherical',bind(BASE/'spherical_sdp/certify.py'))
    nested_module=load('joint_count_fresh_nested',bind(BASE/'spherical_sdp/count_quadratic.py'))
    joint_module=load('joint_count_fresh_joint',bind(BASE/'conditional/joint_count_polynomial.py'))
    affine=load('joint_count_fresh_affine',bind(BASE/'conditional/affine_latitude_counts.py'))
    nested=origin['nested_count_quadratic']
    cone=nested_module.NestedCountCone(nested['radii'],nested['bounds'])
    nested_report=cone.verify(list(map(Q,nested['coefficients'])))
    assert nested['bounds']==[24,25,26]
    assert list(map(Q,nested['radii']))==[breaks[i] for i in [3,4,5]]
    joint_reports=[]
    records=[nested]+origin['conditional_count_polynomials']
    assert len(records)==3
    proof_records=[]
    specifications={}
    for index,record in enumerate(records):
        assert all(Q(r) in breaks for r in record['radii'])
        count_proofs=record['count_proofs']
        assert len(count_proofs)==len(record['radii'])
        for radius,maximum,proof in zip(record['radii'],record['bounds'],count_proofs):
            assert Q(radius)<=Q(proof['R']) and proof['M']<=maximum
            path=bind(BASE/'spherical_sdp'/proof['certificate'])
            assert sha(path)==proof['certificate_sha256']
            proof_records.append(proof)
        if index:
            assert record['bounds']==[24,24,24,25,26]
            assert list(map(Q,record['radii']))==breaks[1:6]
            assert Q(record['weight'])>=0 and Q(record['rounding_origin_repair'])==0
            poly=joint_module.ConditionalCountPolynomial(record)
            assert poly.implications==[(0,4),(1,3)]
            fresh=poly.verify()
            assert fresh['count_vectors_checked']==165503
            joint_reports.append(fresh)
            for implication in record['implications']:
                i,j=implication['inner_index'],implication['outer_index']
                proof=implication['proof']
                key=(record['radii'][i],record['radii'][j],proof['spherical_cutoff'])
                assert implication['maximum']==24 and implication['proof_kind']=='affine_latitude'
                assert Q(proof['inner_radius'])==Q(key[0]) and Q(proof['outer_radius'])==Q(key[1])
                assert proof['maximum_inner_count']==24
                specifications[key]=proof
    checked={}
    for proof in proof_records:
        name=proof['certificate']
        if name not in checked:
            original=BASE/'spherical_sdp'/name
            copy=output_dir/name
            copy.write_bytes(original.read_bytes())
            assert sha(copy)==proof['certificate_sha256']
            checked[name]=spherical.verify(copy)
            assert checked[name]['certificate_sha256']==proof['certificate_sha256']
        report=checked[name]
        r=Q(proof['R'])
        assert report['count_bound']<=proof['M']
        assert report['status']=='EXACT_SPHERICAL_COUNT_'+str(report['count_bound'])+'_VERIFIED'
        assert r>=1 and r*r-r-1<0
        assert 1-Q(1)/(2*r*r)<=Q(report['cosine_cutoff'])
    assert len(checked)==3
    geometry=[]
    for key,proof in specifications.items():
        r,R,s=map(Q,key)
        assert 1<r<R<breaks[-1] and R*R<2
        fresh=affine.check(r,R,s)
        assert fresh['status']=='EXACT_AFFINE_LATITUDE_DOMAIN_CERTIFIED'
        for field in ['inner_cosine_bound','pole_cosine_bound','beta','alpha','transformed_pole_upper']:
            assert Q(fresh[field])==Q(proof[field])
        assert Q(fresh['transformed_pole_upper'])<=s
        dependency=checked[proof['spherical_certificate']]
        assert dependency['certificate_sha256']==proof['spherical_certificate_sha256']
        assert dependency['count_bound']<=24 and s<=Q(dependency['cosine_cutoff'])
        geometry.append(dict(inner_radius=str(r),outer_radius=str(R),proof=fresh,
                             spherical_source_sha256=dependency['certificate_sha256']))
    assert len(geometry)==2
    for path,digest in binding.items():assert sha(ROOT/path)==digest
    return dict(status='FRESH_ACCEPTED_JOINT_COUNT_GEOMETRY_AND_INTEGER_ENERGIES_VERIFIED',
                source_sha256=binding,checker_sha256=sha(Path(__file__)),
                nested_count_polynomial=nested_report,conditional_count_polynomials=joint_reports,
                spherical_counts=checked,affine_implications=geometry,
                proved_count_bounds=[24,24,24,25,26],
                proved_implications=['(23-N0)*(N4-N0)>=0','(23-N1)*(N3-N1)>=0'],
                global_bound_proved=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    report=verify(ROOT/'current/candidates/round9_six_graphs_degree18',
                  ROOT/'audit/round9_six_graphs_monotone_origin_credit_global_verified.json',args.output_dir)
    path=args.output_dir/'FRESH_COUNT_DEPENDENCIES.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],report_sha256=sha(path)),indent=2))
