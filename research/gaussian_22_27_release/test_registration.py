"""Check registration and subprocess wiring; this is not a proof replay."""
from pathlib import Path
import subprocess,sys,json,io,contextlib
R=Path(__file__).resolve().parents[2];orig_run=subprocess.run;orig_read=Path.read_text
out=[]
for n in [22,23,26,27]:
 P=R/f'verifications/dim{n}';script=P/'verify.py'
 for mode in ['valid','bad_bound','bad_package','bad_hash']:
  calls=[]
  def fake(cmd,**kw):calls.append((cmd,str(kw.get('cwd'))));return subprocess.CompletedProcess(cmd,0)
  def read(p,*args,**kw):
   text=orig_read(p,*args,**kw)
   if p==R/'RESULTS.json' and mode in ['bad_bound','bad_package']:
    a=json.loads(text);row=next(x for x in a if x['dimension']==n);row['strict_density_upper_bound' if mode=='bad_bound' else 'package']='0' if mode=='bad_bound' else 'wrong';return json.dumps(a)
   if p==P/'source_hashes.json' and mode=='bad_hash':
    a=json.loads(text);a[next(iter(a))]='0'*64;return json.dumps(a)
   return text
  subprocess.run=fake;Path.read_text=read
  try:
   with contextlib.redirect_stdout(io.StringIO()):exec(compile(orig_read(script),str(script),'exec'),{'__file__':str(script),'__name__':'__main__'})
   assert mode=='valid'
   expected=[([sys.executable,str(P/'verify_primary.py')],str(P))]
   if n in [26,27]:expected+=[([sys.executable,str(P/'independent/verify_all.py')],str(P))]
   else:expected += [([sys.executable,str(P.parent/'dim22_23_independent'/f),str(n)],str(P.parent/'dim22_23_independent'))for f in ['verify_angular.py','verify_radial.py']]
   assert calls==expected
  except AssertionError:
   if mode=='valid':raise
   assert not calls
  finally:subprocess.run=orig_run;Path.read_text=orig_read
 for f in [script,P/'verify_primary.py']:
  result=subprocess.run([sys.executable,'-O',str(f)],capture_output=True,text=True);assert result.returncode!=0
 out.append({'dimension':n,'pass':True,'scope':'Registry/hash rejection and exact subprocess invocation; proof subprocesses mocked.'})
(R/'research/gaussian_22_27_release/REGISTRATION_TESTS.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS four registration gates and negative controls')
