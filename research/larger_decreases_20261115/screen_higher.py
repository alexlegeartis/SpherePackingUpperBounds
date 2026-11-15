"""Numerical scale sweep; every candidate still requires full exact replay."""
from pathlib import Path
from fractions import Fraction as Q
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,os,sys,subprocess,time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ROWS=json.loads((ROOT/'RESULTS.json').read_text())
LOG=HERE/'screen_logs';LOG.mkdir(exist_ok=True)
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')

def run(task):
    n,factor=task
    row=next(r for r in ROWS if r['dimension']==n)
    source=json.loads((ROOT/row['package']/'certificate.json').read_text())
    s=Q(source['s_exact'])*factor
    tag=f'scale_n{n}_f{str(factor).replace("/","_")}'
    cmd=[sys.executable,str(HERE/'search_radial.py'),str(n),'--degree','95',
         '--s',str(s),'--a',source['a_exact'],'--b',source['b_exact'],'--cut','60',
         '--angular',str(ROOT/source['angular_file']),'--tag',tag]
    start=time.time()
    res=subprocess.run(cmd,cwd=ROOT,env=ENV,capture_output=True,text=True)
    (LOG/(tag+'.log')).write_text(res.stdout+res.stderr)
    result=dict(dimension=n,scale_factor=str(factor),returncode=res.returncode,seconds=time.time()-start,tag=tag)
    if res.returncode==0:
        data=json.loads((HERE/tag/'candidate.json').read_text())
        result.update(density_candidate=data['density_candidate'],previous_strict=row['strict_density_upper_bound'],
                      numerical_decrease_percent=100*(1-data['density_candidate']/float(row['strict_density_upper_bound'])))
    return result

if __name__=='__main__':
    tasks=[(n,Q(f)) for n in [29,35,39,45] for f in ['4/5','9/10','1','11/10','6/5']]
    results=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run,t) for t in tasks]):
            result=future.result();results.append(result)
            print(json.dumps(result),flush=True)
            (HERE/'scale_screen.json').write_text(json.dumps(results,indent=2)+'\n')
