"""Registry/hash compatibility and rejection controls, with explicit replay scope."""
import argparse,hashlib,json,runpy,shutil,subprocess,sys,tempfile
from pathlib import Path
from unittest.mock import patch
if not __debug__:raise RuntimeError('Assertions must remain enabled')
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('dimensions',type=int,nargs='+');args=ap.parse_args()
rows=json.loads((ROOT/'RESULTS.json').read_text());records=[]
for n in args.dimensions:
    D=ROOT/f'verifications/dim{n}/refinement_profile20261114'
    cfg=json.loads((D/'config.json').read_text());result=json.loads((D/'full_replay_result.json').read_text())
    row=next(z for z in rows if z['dimension']==n)
    assert row['package']==cfg['package'] and row['strict_density_upper_bound']==cfg['strict_density_upper_bound']
    assert result['pass'] and result['complete_primary_and_independent_replay']
    assert result['strict_density_upper_bound']==cfg['strict_density_upper_bound']
    manifest_sha=hashlib.sha256((D/'source_hashes.json').read_bytes()).hexdigest()
    assert result['source_manifest_sha256']==manifest_sha
    for name,digest in json.loads((D/'source_hashes.json').read_text()).items():
        assert hashlib.sha256((D/name).read_bytes()).hexdigest()==digest,(n,name)
    assert {z['script'] for z in result['executed']}=={'primary_angular.py','primary_radial.py','independent_angular.py','independent_radial.py'}
    for z in result['executed']:assert z['returncode']==0
    for stem in ['primary_angular','primary_radial','independent_angular','independent_radial']:
        gate=json.loads((D/(stem+'_result.json')).read_text());assert gate['pass'] and gate['dimension']==n
    records.append({'dimension':n,'all_four_frozen_proof_paths_pass':True,'registry_matches':True,'source_manifest_sha256':manifest_sha})

changed=subprocess.check_output(['git','diff','--name-only','8a080d6','--','verifications'],cwd=ROOT,text=True).splitlines()
historical=[]
def without_registry_assertions(text):
    return '\n'.join(z for z in text.splitlines() if not(z.startswith('assert ') and "row['" in z))
for name in changed:
    assert name.endswith('/verify.py'),('Historical proof input changed',name)
    old=subprocess.check_output(['git','show','8a080d6:'+name],cwd=ROOT,text=True)
    assert without_registry_assertions(old)==without_registry_assertions((ROOT/name).read_text()),name
    historical.append(ROOT/name)

class ProofCallReached(Exception):pass
def stop_before_proof(*a,**k):raise ProofCallReached
wrappers=historical+[ROOT/f'verifications/dim{n}/refinement_profile20261114/verify.py' for n in args.dimensions]
compatibility=[]
for file in wrappers:
    reached=False
    with patch('subprocess.run',side_effect=stop_before_proof):
        try:runpy.run_path(str(file))
        except ProofCallReached:reached=True
    assert reached,file
    compatibility.append({'wrapper':str(file.relative_to(ROOT)),'registry_and_hash_prefix_pass':True,'scope':'Stopped immediately before first proof subprocess; not a new historical mathematical replay'})

controls=[]
def rejection(command,label):
    result=subprocess.run(command,capture_output=True,text=True,cwd=ROOT)
    assert result.returncode!=0,(label,result.stdout,result.stderr)
    assert 'AssertionError' in result.stderr or 'Assertions' in result.stderr,(label,result.stderr)
    controls.append({'control':label,'rejected':True,'returncode':result.returncode})
for file in wrappers:
    rejection([sys.executable,'-O',str(file)],'optimized assertions: '+str(file.relative_to(ROOT)))
shared=ROOT/'verifications/profile_enrichment_20261114';original=ROOT/'verifications/dim38/refinement_profile20261114'
for stem in ['primary_angular','primary_radial','independent_angular','independent_radial']:
    rejection([sys.executable,'-O',str(shared/(stem+'.py')),str(original)],'optimized assertions: '+stem)
with tempfile.TemporaryDirectory(prefix='negative_profile_',dir=original.parent) as temporary:
    D=Path(temporary)
    for name in ['certificate.json','angular.json','config.json','verify.py','source_hashes.json']:
        shutil.copy2(original/name,D/name)
    cert=json.loads((D/'certificate.json').read_text());cert['fourier_polynomial'][0]='12345'
    (D/'certificate.json').write_text(json.dumps(cert))
    rejection([sys.executable,str(D/'verify.py')],'tampered certificate rejected by manifest')
    for stem in ['primary_radial','independent_radial']:
        rejection([sys.executable,str(shared/(stem+'.py')),str(D)],'wrong Fourier coefficient: '+stem)
    shutil.copy2(original/'certificate.json',D/'certificate.json')
    angular=json.loads((D/'angular.json').read_text());angular['gegenbauer_coefficients'][1]='-1'
    (D/'angular.json').write_text(json.dumps(angular))
    for stem in ['primary_angular','independent_angular']:
        rejection([sys.executable,str(shared/(stem+'.py')),str(D)],'negative harmonic coefficient: '+stem)
    cert=json.loads((D/'certificate.json').read_text());cert['a_exact']='-121'
    angular['profile_a_exact']='-121'
    (D/'certificate.json').write_text(json.dumps(cert));(D/'angular.json').write_text(json.dumps(angular))
    for stem in ['primary_radial','independent_radial']:
        rejection([sys.executable,str(shared/(stem+'.py')),str(D)],'inadmissible logarithmic slope: '+stem)

# Reject an unrelated registry path even when its bound is unchanged.
read_text=Path.read_text
def bad_registry(self,*a,**k):
    text=read_text(self,*a,**k)
    if self.resolve()==(ROOT/'RESULTS.json').resolve():
        data=json.loads(text)
        next(z for z in data if z['dimension']==38)['package']='verifications/unrelated'
        return json.dumps(data)
    return text
with patch.object(Path,'read_text',bad_registry),patch('subprocess.run',side_effect=stop_before_proof):
    rejected=False
    try:runpy.run_path(str(original/'verify.py'))
    except AssertionError:rejected=True
    assert rejected
controls.append({'control':'unrelated registered package','rejected':True})
rejected_real=HERE/'dim44_rejected_angular'
rejection([sys.executable,str(shared/'primary_angular.py'),str(rejected_real)],'retained genuine 44D angular sign island')
subprocess.run([sys.executable,str(ROOT/'comparison.py')],cwd=ROOT,check=True)
out={'pass':True,'base_commit':'8a080d6a00909391a91ecbac097e4004c404c40d','new_packages':records,
     'historical_changes_only_registry_assertions':True,'historical_mathematical_inputs_and_checks_unchanged':True,
     'wrapper_compatibility':compatibility,'negative_controls':controls,
     'generated_comparisons_and_70_100_digit_agreement':True,
     'scope':'All new mathematical proof paths replayed; changed historical wrapper registry/hash prefixes checked. The full historical run_all.py aggregate was not rerun in this campaign.'}
(HERE/'CAMPAIGN_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
print('CAMPAIGN AUDIT PASS',len(records),'new packages,',len(compatibility),'wrapper prefixes,',len(controls),'rejection controls',flush=True)
