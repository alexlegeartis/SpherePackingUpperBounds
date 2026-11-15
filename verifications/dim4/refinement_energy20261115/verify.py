#!/usr/bin/env python3
"""Replay the full pinned D4 proof and the new Gram-assisted energy certificate.

Default: fresh full base replay, followed by the extension. --base-replay DIR
reuses a separately completed full run and records that reuse explicitly.
--extension-only checks the extension against committed base data; it is not
a fresh replay of the underlying global certificate.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, tempfile, time

HERE=Path(__file__).resolve().parent
VENDOR=HERE/'vendor'
D4=VENDOR/'D4'
UPSTREAM_COMMIT='bb97550ba47ad92bacfee567843888be16895a23'
OLD_BOUND=Q('0.629521681653')
NEW_BOUND=Q('0.629520189662')
sys.dont_write_bytecode=True


def read(path):
    def unique(pairs):
        result={}
        for k,v in pairs:
            if k in result: raise ValueError('duplicate JSON key: '+k)
            result[k]=v
        return result
    return json.loads(Path(path).read_text(),object_pairs_hook=unique)


def exact(value):
    if type(value) not in (str,int): raise ValueError('inexact rational input')
    return Q(value)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_upstream():
    manifest=read(VENDOR/'UPSTREAM.json')
    assert manifest['commit']==UPSTREAM_COMMIT
    assert manifest['repository']=='https://github.com/Felpix-Studios/sphere-packing'
    assert manifest['license']=='Apache-2.0'
    seen=set()
    for row in manifest['files']:
        relative=Path(row['path'])
        assert not relative.is_absolute() and '..' not in relative.parts
        assert row['path'] not in seen
        seen.add(row['path'])
        p=VENDOR/relative
        assert p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],str(p)
    assert len(seen)==278
    return len(seen)


def check_base(base):
    summary=read(base/'SUMMARY.json')
    assert summary['status']=='PASS' and summary['mode']=='full'
    assert summary['full_certificate_replay'] is True and summary['branches_rerun']==112
    stages=['stage1_v7','stage2_origin','stage3_census','stage4_graph_slack',
            'stage5_kernelbox','kernel_independent','kernel_link']
    for stage in stages:
        r=read(base/stage/'RESULT.json')
        assert r['status']=='PASS' and r['stage']==stage
    stage5=read(base/'stage5_kernelbox/RESULT.json')
    assert exact(stage5['lam5'])==Q(1,10**6)
    assert stage5['kernel_box']['j5_roots']==2350 and stage5['kernel_box']['j5_leaves_replayed']>0
    assert exact(summary['graph_saving'])==Q(3,8000000)
    assert exact(summary['U_raw'])==exact(read(D4/'data/census_floors.json')['exact_raw_objective'])
    assert exact(summary['U_new'])<OLD_BOUND
    assert summary['group_floors']==124
    return summary


def check_extension(census_file,workers):
    # Import only after the caller has selected the fresh stage-2 input.
    import polynomial
    import gram
    import boxes
    registry=read(HERE.parents[2]/'RESULTS.json')
    registered=[row for row in registry if row['dimension']==4]
    assert len(registered)==1
    assert exact(registered[0]['strict_density_upper_bound'])==NEW_BOUND
    assert exact(registered[0]['previous_upper_density'])==OLD_BOUND
    assert registered[0]['package']=='verifications/dim4/refinement_energy20261115'
    assert registered[0]['source_commit']==UPSTREAM_COMMIT
    data=read(HERE/'polynomials.json')
    a=exact(data['a'])
    assert a==Q(133,250) and 2+2*a>Q(49,16)
    p1=list(map(exact,data['1']['coefficients']))
    p5=list(map(exact,data['5']['coefficients']))
    assert len(p1)==17 and data['5']['inherited'] is True
    expected5=[Q(1,10**6)*(Q(1,8)-a*a/4),Q(0),
               Q(1,10**6)*3*(Q(3,16)-a*a/4),Q(0),Q(5,16*10**6)]
    assert p5==expected5
    r1=Q(501473,500000);r5=Q(102597,100000)
    assert boxes.kf.BREAKS[1]==r1 and boxes.kf.BREAKS[5]==r5
    central={}
    for j,r,p in [(1,r1,p1),(5,r5,p5)]:
        cap=1-Q(1)/(2*r*r)
        assert 1<=r<Q(3,2) and r/2<=cap<=a
        central[str(j)]=polynomial.check_nonpositive(p,a,cap)
    credits=gram.check_credits(D4/'data/v7/kernel.json',HERE/'gram_credits.json',1,p1)
    assert set(credits['credits'])=={5,6,7,8,9}
    box=boxes.check_tree(HERE/'kernel_box_j1.json',1,p1,workers)
    census=read(census_file)
    frozen=read(D4/'data/census_floors.json')
    assert census['profile']==frozen['profile']
    for k in ['exact_raw_objective','graph_saving','new_beta']:
        assert exact(census[k])==exact(frozen[k])
    expected={(0,1)}|{(n,j) for n in range(1,25) for j in range(1,6)}|{(25,4),(25,5),(26,5)}
    rows={}
    for group in census['profile']:
        n,j=group['N'],group['first_radius']
        assert type(n) is int and type(j) is int and (n,j) not in rows
        floor=exact(group['new_floor']);assert floor>=0
        gain=max(polynomial.energy(p1 if j==1 else p5,n),Q(0))
        rows[n,j]=(floor+gain,gain)
    assert set(rows)==expected
    binding=min(rows,key=lambda key:rows[key][0])
    assert binding==(24,1)
    floor=rows[binding][0]
    raw=exact(census['exact_raw_objective'])
    saving=exact(census['graph_saving']);assert saving==Q(3,8000000)
    upper=raw-floor/8-saving
    assert upper<NEW_BOUND<OLD_BOUND
    grid=Q(-((-upper.numerator*10**12)//upper.denominator),10**12)
    assert grid==NEW_BOUND and 0<grid-upper<Q(1,10**12)
    frozen_result=read(HERE/'certificate.json')
    assert exact(frozen_result['exact_density_upper_bound'])==upper
    assert exact(frozen_result['strict_density_upper_bound'])==NEW_BOUND
    assert exact(frozen_result['binding_energy_gain'])==rows[binding][1]
    for name in ['polynomials.json','gram_credits.json','kernel_box_j1.json','vendor/UPSTREAM.json']:
        assert frozen_result['input_sha256'][name]==sha(HERE/name)
    return dict(strict_density_upper_bound='0.629520189662',exact_density_upper_bound=str(upper),
                old_strict_density_upper_bound='0.629521681653',
                decrease=str(OLD_BOUND-NEW_BOUND),binding_group=list(binding),
                energy_gain=str(rows[binding][1]),
                gram_range_equations_checked=credits['range_equations_checked'],
                central_sign=central,kernel_box_j1=box,groups_checked=len(rows))


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--workers',type=int,choices=range(1,5),default=2)
    scope=parser.add_mutually_exclusive_group()
    scope.add_argument('--base-replay',type=Path)
    scope.add_argument('--extension-only',action='store_true')
    args=parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):parser.error('assertions must be enabled')
    if sys.version_info<(3,12):parser.error('Python 3.12 or later is required')
    if args.out:
        out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    else:out=Path(tempfile.mkdtemp(prefix='packing-d4-energy-'))
    start=time.time()
    result=dict(status='FAIL',full_base_replay=False,base_replay_reused=bool(args.base_replay),
                extension_only=args.extension_only)
    try:
        result['upstream_files_checked']=check_upstream()
        base=args.base_replay.resolve() if args.base_replay else out/'base'
        if not (args.base_replay or args.extension_only):
            env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            env.setdefault('CXX','g++')
            with (out/'base.log').open('w') as log:
                subprocess.run([sys.executable,str(D4/'verify/verify_all.py'),'--out',str(base),
                                '--mode','full','--workers',str(args.workers)],cwd=D4,env=env,
                               stdout=log,stderr=subprocess.STDOUT,check=True)
            result['full_base_replay']=True
        census=D4/'data/census_floors.json'
        if not args.extension_only:
            summary=check_base(base)
            result['base_summary_sha256']=sha(base/'SUMMARY.json')
            result['base_branches_replayed']=summary['branches_rerun']
            os.environ['D4_STAGE2_DIR']=str(base/'stage2_origin')
            census=base/'stage3_census/FINITE_COUNT_ACCEPTANCE.json'
        os.environ['D4_CENSUS_FILE']=str(census)
        result.update(check_extension(census,args.workers))
        result['status']='PASS_EXTENSION_ONLY' if args.extension_only else 'PASS'
    except Exception as exc:
        result['error']=repr(exc)
        raise
    finally:
        result['seconds']=round(time.time()-start,3)
        (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
