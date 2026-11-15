"""Finite exact-root/pair-guided continuation of a retained radial candidate."""
import argparse,json,subprocess,sys,shutil,time
from pathlib import Path
from common import HERE,ROOT
ap=argparse.ArgumentParser();ap.add_argument('candidate');ap.add_argument('--rounds',type=int,default=14)
args=ap.parse_args();file=Path(args.candidate).resolve();D=file.parent
seed=json.loads(file.read_text());base=D/'verification_history';base.mkdir(exist_ok=True)
run=0
while (base/f'run{run:03}').exists():run+=1
H=base/f'run{run:03}';H.mkdir()
history=[];start=time.time()
for i in range(args.rounds):
    R=H/f'round{i:02}';R.mkdir()
    cutsfile=D/'cuts.json';cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {'physical':[],'fourier':[],'pair':[]}
    if (D/'root_rejection.json').exists():
        for kind,roots in json.loads((D/'root_rejection.json').read_text())['global_sign_roots'].items():
            xs=sorted(z[0] for z in roots);cuts[kind]+=xs+[(x+y)/2 for x,y in zip(xs,xs[1:])]
    if (D/'pair_rejection.json').exists():cuts['pair']+=json.loads((D/'pair_rejection.json').read_text())['pair_points']
    cuts={k:sorted(set(v)) for k,v in cuts.items()};cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
    cmd=[sys.executable,str(HERE/'search_radial.py'),str(seed['dimension']),'--degree',str(seed['degree']),
         '--s',seed['s_exact'],'--a',seed['a_exact'],'--b',seed['b_exact'],'--q',seed['q_exact'],
         '--cut',seed['cut_X_exact'],'--angular',str(ROOT/seed['angular_file']),'--tag',D.name]
    search=subprocess.run(cmd,capture_output=True,text=True)
    (R/'search.log').write_text(search.stdout+search.stderr);print(search.stdout,flush=True)
    if search.returncode:print(search.stderr,flush=True);raise RuntimeError('Search failed')
    shutil.copy2(D/'candidate.json',R/'candidate.json');shutil.copy2(cutsfile,R/'cuts.json')
    check=subprocess.run([sys.executable,str(HERE/'check_radial.py'),str(D/'candidate.json')],capture_output=True,text=True)
    (R/'check.log').write_text(check.stdout+check.stderr);print(check.stdout,flush=True)
    for name in ['radial_export.json','root_rejection.json','pair_rejection.json','radial_certificate.json']:
        if (D/name).exists():shutil.copy2(D/name,R/name)
    row={'round':i,'returncode':check.returncode,'density_candidate':json.loads((D/'candidate.json').read_text())['density_candidate'],'seconds':time.time()-start}
    history.append(row);(H/'history.json').write_text(json.dumps(history,indent=2)+'\n');print(row,flush=True)
    if check.returncode==0:break
    if not any(z in check.stderr for z in ['Global radial sign gate failed','Pair gate failed','Pair interval proof inconclusive']):
        print(check.stderr,flush=True);raise RuntimeError('Unexpected exact-check failure')
else:raise RuntimeError('Finite continuation budget exhausted')
