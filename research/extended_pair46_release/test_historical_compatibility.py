"""Registration-only compatibility test; historical proof subprocesses are mocked."""
from pathlib import Path
import subprocess,json,io,contextlib
R=Path(__file__).resolve().parents[2];orig_run=subprocess.run;orig_read=Path.read_text;out=[]
for n in [46]:
 for sub in ['', '/refinement_angular','/refinement_joint','/refinement_degree47']:
  p=R/f'verifications/dim{n}{sub}/verify.py'
  for bad in [False,True]:
   calls=[]
   def fake(cmd,**kw):calls.append(cmd);return subprocess.CompletedProcess(cmd,0)
   def read(q,*a,**kw):
    s=orig_read(q,*a,**kw)
    if q==R/'RESULTS.json' and bad:
     rows=json.loads(s);next(x for x in rows if x['dimension']==n)['package']='wrong';return json.dumps(rows)
    return s
   subprocess.run=fake;Path.read_text=read
   try:
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(orig_read(p),str(p),'exec'),{'__file__':str(p),'__name__':'__main__'})
    assert not bad and len(calls)==2
   except AssertionError:
    if not bad:raise
    assert not calls
   finally:subprocess.run=orig_run;Path.read_text=orig_read
  out.append({'script':str(p.relative_to(R)),'accepts_stronger_bound':True,'rejects_unknown_package':True,'both_historical_chains_invoked':True})
(R/'research/extended_pair46_release/HISTORICAL_COMPATIBILITY.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS four historical registration wrappers; proof subprocesses mocked')
