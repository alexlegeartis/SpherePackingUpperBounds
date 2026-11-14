"""Collect exact active contact rays with randomized reduced bases (capped)."""
import json
import pathlib
import subprocess
import sys
import numpy as np

root=pathlib.Path(__file__).resolve().parent;p=root/'dim48mod3'
lo=int(sys.argv[1]);hi=int(sys.argv[2]);cap=int(sys.argv[3]) if len(sys.argv)>3 else 300000000
reports=[]
for seed in range(lo,hi):
    subprocess.run([sys.executable,str(root/'periodic44_setup.py'),str(seed)],check=True)
    output=p/f'periodic44_raw{seed}.txt'
    result=subprocess.run([str(root.parent/'bin'/'periodic44_shell'),'48',str(p/'periodic44_Q.txt'),str(p/'periodic44_R.txt'),'6',str(cap),str(output),'0',str(p/'periodic44_w.txt')],capture_output=True,text=True,check=True)
    R=np.loadtxt(output,dtype=np.int64,ndmin=2);T=np.load(p/f'periodic44_basis{seed}.npz')['T']
    np.save(p/f'periodic44_shell{seed}.npy',R[:,1:]@T)
    report={'seed':seed,**json.loads(result.stdout),'integer_vectors_saved':len(R)}
    reports.append(report);print(json.dumps(report),flush=True)
    (p/f'periodic44_sampling_{lo}_{hi}.json').write_text(json.dumps(reports,indent=2))
