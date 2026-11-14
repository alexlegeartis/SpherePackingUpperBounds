"""Register only packages with successful complete frozen replays."""
import argparse,hashlib,json,subprocess
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ap=argparse.ArgumentParser();ap.add_argument('dimensions',type=int,nargs='+');args=ap.parse_args()
rows=json.loads((ROOT/'RESULTS.json').read_text())
baseline=json.loads(subprocess.check_output(['git','show','8a080d6:RESULTS.json'],cwd=ROOT,text=True))
for n in args.dimensions:
    D=ROOT/f'verifications/dim{n}/refinement_profile20261114'
    cfg=json.loads((D/'config.json').read_text());result=json.loads((D/'full_replay_result.json').read_text())
    assert result['pass'] and result['complete_primary_and_independent_replay'] and result['dimension']==n
    assert result['strict_density_upper_bound']==cfg['strict_density_upper_bound']
    assert len(result['executed'])==4 and all(z['returncode']==0 for z in result['executed'])
    assert result['source_manifest_sha256']==hashlib.sha256((D/'source_hashes.json').read_bytes()).hexdigest()
    for name,digest in json.loads((D/'source_hashes.json').read_text()).items():
        assert hashlib.sha256((D/name).read_bytes()).hexdigest()==digest,name
    row=next(z for z in rows if z['dimension']==n);old=next(z for z in baseline if z['dimension']==n)
    assert Fraction(cfg['strict_density_upper_bound'])<Fraction(old['strict_density_upper_bound'])
    previous=row.get('earlier_repository_upper_densities',[])
    if old['strict_density_upper_bound'] not in previous:previous.append(old['strict_density_upper_bound'])
    packages=row.get('earlier_repository_packages',[])
    candidates=[old['package'],f'verifications/dim{n}']
    candidates+=[str(p.parent.relative_to(ROOT)) for p in (ROOT/f'verifications/dim{n}').glob('refinement_*/verify.py')]
    for name in candidates:
        if name!=cfg['package'] and name not in packages:packages.append(name)
    row.update(strict_density_upper_bound=cfg['strict_density_upper_bound'],package=cfg['package'],
               source_package=cfg['package'],previous_repository_package=old['package'],
               previous_repository_strict_upper=old['strict_density_upper_bound'],
               earlier_repository_upper_densities=previous,earlier_repository_packages=packages,
               method='Signed logarithmic radial profile, new dimension-specific angular kernel and degree95 Laguerre enrichment; complete primary and independent continuum proofs',
               radial_profile={'a_exact':cfg['profile_a_exact'],'b_exact':cfg['profile_b_exact'],'q_exact':'1'},
               external_neighbor_counts_required=False)
    wrappers=[ROOT/f'verifications/dim{n}/verify.py']
    wrappers+=[p for p in (ROOT/f'verifications/dim{n}').glob('refinement_*/verify.py') if p.parent!=D]
    for file in wrappers:
        text=file.read_text();lines=text.splitlines(keepends=True)
        for i,line in enumerate(lines):
            if "row['package']" not in line or not line.lstrip().startswith('assert '):continue
            if "==cfg['package']" in line:
                line=line.replace("row['package']==cfg['package']",f"row['package'] in [cfg['package'], '{cfg['package']}']")
            elif " in [" in line:
                if repr(cfg['package']) not in line:line=line.rstrip('\n')[:-1]+f", '{cfg['package']}']\n"
            else:
                needle=f"row['package']=='verifications/dim{n}'"
                assert needle in line,(file,line)
                line=line.replace(needle,f"row['package'] in ['verifications/dim{n}', '{cfg['package']}']")
            line=line.replace("Fraction(row['strict_density_upper_bound'])==", "Fraction(row['strict_density_upper_bound'])<=")
            lines[i]=line
        # The earlier radial wrappers put package and bound checks on separate lines.
        for i,line in enumerate(lines):
            if line.startswith("assert Fraction(row['strict_density_upper_bound'])==Fraction(cfg["):
                lines[i]=line.replace(')==Fraction(',')<=Fraction(')
        file.write_text(''.join(lines))
    print(n,old['strict_density_upper_bound'],'->',cfg['strict_density_upper_bound'])
(ROOT/'RESULTS.json').write_text(json.dumps(rows,indent=2)+'\n')
runall=ROOT/'run_all.py';text=runall.read_text()
line="jobs += [(f'dim{n}', 'refinement_profile20261114/verify.py') for n in "+repr(sorted(args.dimensions))+']\n'
assert 'refinement_profile20261114/verify.py' not in text,'Run this registrar once for the complete selected set'
text=text.replace("jobs.append(('dim7','verify.py'))",line+"jobs.append(('dim7','verify.py'))")
runall.write_text(text)
comparison=ROOT/'comparison.py';text=comparison.read_text()
needle='assert set(proof_ideas)'
position=text.index(needle)
text=text[:position]+"for n in "+repr(sorted(args.dimensions))+":proof_ideas[n]='Signed logarithmic profile + degree95 radial enrichment'\n"+text[position:]
comparison.write_text(text)
