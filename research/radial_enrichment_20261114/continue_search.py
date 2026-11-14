"""Retain finite search rounds; root isolation and interval proof decide acceptance."""
import argparse, json, subprocess, sys, shutil, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb_poly, arb, fmpq as Q, ctx

ap=argparse.ArgumentParser()
ap.add_argument('dimension',type=int)
ap.add_argument('--degree',type=int,default=47)
ap.add_argument('--cut',type=float,default=40)
ap.add_argument('--rounds',type=int,default=16)
ap.add_argument('--margin',type=float,default=1e-7)
args=ap.parse_args();ctx.prec=512
HERE=Path(__file__).resolve().parent
D=HERE/f'n{args.dimension}_d{args.degree}_c{args.cut:g}'
archive=D/'history';archive.mkdir(exist_ok=True)
run=0
while (archive/f'run{run:03}').exists():run+=1
H=archive/f'run{run:03}';H.mkdir()
history=[];start=time.time()

def add_cuts():
    file=D/'cuts.json'
    cuts=json.loads(file.read_text()) if file.exists() else {'physical':[],'fourier':[],'pair':[]}
    if (D/'root_result.json').exists():
        bad=json.loads((D/'root_result.json').read_text())['global_sign_roots']
        for kind,roots in bad.items():
            xs=sorted(v[0] for v in roots)
            cuts[kind]+=xs+[(a+b)/2 for a,b in zip(xs,xs[1:])]
    if (D/'export.json').exists():
        data=json.loads((D/'export.json').read_text())
        poly=arb_poly([arb(Q(c)) for c in data['physical_polynomial']])
        M=arb(Q(data['M_exact']));pi=arb.pi()
        def value(x):
            xx=arb(float(x));hv=(-Q(17,8)*(xx-6*pi)).exp()
            return float((poly(xx)*(Q(9,8)*xx-Q(51,4)*pi).exp()/(M*(2-hv))).mid())
        grid=np.linspace(6*np.pi,args.cut,3001)
        vals=np.array([value(x) for x in grid])
        peaks=[]
        for i in range(1,len(grid)-1):
            if vals[i]>vals[i-1] and vals[i]>vals[i+1]:
                result=minimize_scalar(lambda x:-value(x),bounds=(grid[i-1],grid[i+1]),
                    method='bounded',options={'xatol':1e-12})
                if -result.fun>0.999:
                    peaks.append((float(result.x),float(-result.fun)))
        cuts['pair'] += [x for x,v in peaks]
        print('PAIR_PEAKS',peaks,flush=True)
    cuts={k:sorted(set(v)) for k,v in cuts.items()}
    file.write_text(json.dumps(cuts,indent=2)+'\n')

for i in range(args.rounds):
    add_cuts()
    rundir=H/f'round{i:02}';rundir.mkdir(exist_ok=True)
    command=[sys.executable,str(HERE/'search.py'),str(args.dimension),'--degree',str(args.degree),
             '--cut',str(args.cut),'--margin',str(args.margin)]
    result=subprocess.run(command,capture_output=True,text=True)
    (rundir/'search.log').write_text(result.stdout+result.stderr)
    print(result.stdout,flush=True)
    if result.returncode:
        print(result.stderr,flush=True);raise RuntimeError('search failed')
    shutil.copy2(D/'candidate.json',rundir/'candidate.json')
    shutil.copy2(D/'cuts.json',rundir/'cuts.json')
    audit=subprocess.run([sys.executable,str(HERE/'export_check.py'),str(D/'candidate.json')],capture_output=True,text=True)
    (rundir/'audit.log').write_text(audit.stdout+audit.stderr)
    print(audit.stdout,flush=True)
    for name in ['root_result.json','export.json','pair_failure.json','radial_certificate.json']:
        if (D/name).exists():shutil.copy2(D/name,rundir/name)
    row={'round':i,'density_candidate':json.loads((D/'candidate.json').read_text())['density_candidate'],
         'radial_audit_returncode':audit.returncode,'seconds':time.time()-start}
    if audit.returncode and not any(x in audit.stderr for x in ['Forbidden sign roots','Pair coupling inconclusive']):
        print(audit.stderr,flush=True);raise RuntimeError('unexpected audit failure')
    history.append(row)
    (H/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    (D/'history.json').write_text(json.dumps(history,indent=2)+'\n')
    print(row,flush=True)
    if audit.returncode==0:
        print('ACCEPTED RADIAL GATES; ANGULAR REPLAY STILL REQUIRED',flush=True)
        break
else:
    raise RuntimeError('finite repair budget exhausted')
