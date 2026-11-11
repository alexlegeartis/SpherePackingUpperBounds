from pathlib import Path
import subprocess,threading,time,sys,contextlib,io,json
from collections import Counter
R=Path(__file__).resolve().parents[2];original=subprocess.run

def exercise(workers,fail=False):
 calls=[];active={};lock=threading.Lock();state={'current':0,'maximum':0};out=io.StringIO()
 def fake(cmd,**kw):
  cwd=str(kw.get('cwd',''));script=cmd[-1]
  with lock:
   if cwd:assert active.get(cwd,0)==0;active[cwd]=1
   state['current']+=1;state['maximum']=max(state['maximum'],state['current']);calls.append((tuple(cmd),cwd))
  time.sleep(.002)
  with lock:
   if cwd:active[cwd]=0
   state['current']-=1
  if fail and cwd.endswith('dim25'):raise subprocess.CalledProcessError(1,cmd)
  return subprocess.CompletedProcess(cmd,0)
 subprocess.run=fake;oldargv=sys.argv;sys.argv=[str(R/'run_all.py'),'--jobs',str(workers)]
 try:
  with contextlib.redirect_stdout(out):exec(compile((R/'run_all.py').read_text(),str(R/'run_all.py'),'exec'),{'__file__':str(R/'run_all.py'),'__name__':'__main__'})
  assert not fail
 except subprocess.CalledProcessError:assert fail
 finally:subprocess.run=original;sys.argv=oldargv
 if fail:
  assert not any('origin_slack/verify_algebra.py'in c[0][-1]for c in calls)
  assert 'PASS: all registered' not in out.getvalue()
 else:assert 'PASS: all registered proof replays'in out.getvalue()
 return calls,state
serial,s=exercise(1);parallel,p=exercise(4);assert Counter(serial)==Counter(parallel);assert p['maximum']>1
exercise(4,True)
out={'pass':True,'same_full_command_inventory':True,'serial_call_count':len(serial),'parallel_peak':p['maximum'],'same_dimension_overlap':False,'failure_prevents_final_shared_checks':True,'note':'Runner behavior test with mocked subprocesses; not a mathematical proof replay.'}
(R/'research/gaussian21_release/RUNNER_TESTS.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
