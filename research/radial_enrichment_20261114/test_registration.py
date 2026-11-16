"""Registry, immutable inputs, complete subprocess coverage and -O rejection."""
import contextlib,io,json,subprocess,sys
from pathlib import Path
from fractions import Fraction

ROOT=Path(__file__).resolve().parents[2]
original_run=subprocess.run;original_read=Path.read_text
out=[]
for n in [45,46,47]:
    D=ROOT/f'verifications/dim{n}/refinement_radial20261114'
    script=D/'verify.py'
    frozen_bound=Fraction(json.loads(original_read(D/'config.json'))['strict_density_upper_bound'])
    for mode in ['valid','bad_bound','bad_package','bad_hash']:
        calls=[]
        def fake(cmd,**kwargs):
            calls.append(cmd[-1]);return subprocess.CompletedProcess(cmd,0,stdout='',stderr='')
        def read(path,*args,**kwargs):
            s=original_read(path,*args,**kwargs)
            if path==ROOT/'RESULTS.json' and mode in ['bad_bound','bad_package']:
                rows=json.loads(s);r=next(r for r in rows if r['dimension']==n)
                # Historical certificates accept an already registered stronger
                # successor. A smaller bound is therefore a valid compatibility
                # input; use a weaker bound to test the rejection gate.
                r['strict_density_upper_bound' if mode=='bad_bound' else 'package']=str(2*frozen_bound) if mode=='bad_bound' else 'wrong'
                return json.dumps(rows)
            if path==D/'source_hashes.json' and mode=='bad_hash':
                r=json.loads(s);r[next(iter(r))]='0'*64;return json.dumps(r)
            return s
        subprocess.run=fake;Path.read_text=read
        # Keep generated results untouched: mocked subprocesses only test the gates.
        original_write=Path.write_text
        Path.write_text=lambda *a,**kw:0
        rejected=False
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(original_read(script),str(script),'exec'),{'__file__':str(script),'__name__':'__main__'})
        except AssertionError:rejected=True
        finally:
            subprocess.run=original_run;Path.read_text=original_read;Path.write_text=original_write
        assert rejected==(mode!='valid')
        if mode=='valid':
            assert [Path(p).name for p in calls]==[
                'independent_cyclic40_angular.py','independent_cyclic40_domain.py',
                'verify_radial_primary.py','audit_angular_cyclic.py','audit_radial.py']
        else:assert not calls
    for p in [script,D/'verify_radial_primary.py',D/'independent/audit_angular_cyclic.py',D/'independent/audit_radial.py']:
        r=original_run([sys.executable,'-O',str(p)],capture_output=True,text=True)
        assert r.returncode!=0
    # Compatibility gates of every preserved 45/46/47 proof chain.
    paths=[ROOT/f'verifications/dim{n}/verify.py']
    if n in [46,47]:
        paths += [ROOT/f'verifications/dim{n}'/p/'verify.py' for p in ['refinement_angular','refinement_joint','refinement_degree47']]
    if n==46:paths.append(ROOT/'verifications/dim46/refinement_extended_pair/verify.py')
    for p in paths:
        calls=[];subprocess.run=fake
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(original_read(p),str(p),'exec'),{'__file__':str(p),'__name__':'__main__'})
        finally:subprocess.run=original_run
        assert len(calls)==2,(p,calls)
    out.append({'dimension':n,'valid_registry_runs_all_five_proofs':True,
                'weaker_bound_bad_package_and_bad_hash_rejected_before_subprocess':True,
                'all_new_entrypoints_reject_optimized_mode':True,
                'preserved_historical_registration_gates_pass':True,
                'scope':'Control-flow gates; mathematical replays are recorded separately'})
(Path(__file__).parent/'REGISTRATION_TESTS.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS registry, hashes, complete gate coverage, historical compatibility and -O rejection')
