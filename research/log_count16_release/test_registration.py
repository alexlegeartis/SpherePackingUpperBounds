from pathlib import Path
import subprocess,sys,json,io,contextlib,hashlib
R=Path(__file__).resolve().parents[2];orig_run=subprocess.run;orig_read=Path.read_text;rows=json.loads((R/'RESULTS.json').read_text());out=[]
for n in [16]:
 P=R/f'verifications/dim{n}';script=P/'verify.py'
 for mode in ['valid','bad_bound','bad_package','bad_hash']:
  calls=[]
  def fake(cmd,**kw):calls.append((cmd,str(kw.get('cwd'))));return subprocess.CompletedProcess(cmd,0)
  def read(p,*args,**kw):
   text=orig_read(p,*args,**kw)
   if p==R/'RESULTS.json' and mode in ['bad_bound','bad_package']:
    a=json.loads(text);r=next(x for x in a if x['dimension']==n);r['strict_density_upper_bound' if mode=='bad_bound' else 'package']='0' if mode=='bad_bound' else 'wrong';return json.dumps(a)
   if p==P/'source_hashes.json'and mode=='bad_hash':
    a=json.loads(text);a[next(iter(a))]='0'*64;return json.dumps(a)
   return text
  subprocess.run=fake;Path.read_text=read
  try:
   with contextlib.redirect_stdout(io.StringIO()):exec(compile(orig_read(script),str(script),'exec'),{'__file__':str(script),'__name__':'__main__'})
   assert mode=='valid'
   assert calls==[([sys.executable,str(P/f)],str(P))for f in ['verify_all.py','independent/dim16_log_count_combined_independent_20261111/run_all.py']]
  except AssertionError:
   if mode=='valid':raise
   assert not calls
  finally:subprocess.run=orig_run;Path.read_text=orig_read
 for script2 in [script,P/'verify_all.py',P/'independent/dim16_log_count_combined_independent_20261111/run_all.py']:
  c=subprocess.run([sys.executable,'-O',str(script2)],capture_output=True,text=True);assert c.returncode!=0
 out.append({'dimension':n,'registry_and_sourcehash_gates':'PASS','both_complete_chains_invoked':'PASS','incorrect_bound_package_and_hash_rejected':'PASS','optimized_mode_rejected':'PASS','scope':'Changed registration behavior tested with mocked proof subprocesses; completed mathematical audits are separately hash-bound.'})
(R/'research/log_count16_release/REGISTRATION_TESTS.json').write_text(json.dumps(out,indent=2)+'\n');print('11–14D registration, negative-binding and optimized-mode tests PASS')
