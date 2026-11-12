"""Verify rejection of deliberately corrupted copies; original inputs unchanged."""
from pathlib import Path
import tempfile,shutil,json,subprocess,sys
D=Path(__file__).parent;results={}
for name in ['fourier_mismatch','negative_g','angular_negative','component_mismatch','count_bound','negative_count_coefficient','optimized_execution']:
 with tempfile.TemporaryDirectory() as t:
  p=Path(t)/'package';shutil.copytree(D,p,ignore=shutil.ignore_patterns('__pycache__'))
  if name=='optimized_execution':
   cmd=[sys.executable,'-O',str(p/'verify_all.py')]
  else:
   if name=='fourier_mismatch':file=p/'radial/certificate.json';j=json.loads(file.read_text());j['fourier'][1]='0';script='radial/verify.py'
   elif name=='negative_g':file=p/'radial/fit.json';j=json.loads(file.read_text());j['parameters'][0]=-1000000000;script='radial/verify.py'
   elif name=='angular_negative':file=p/'geometry/angular.json';j=json.loads(file.read_text());j['gegenbauer_coefficients'][0]='-1';script='geometry/verify_structure.py'
   elif name=='count_bound':file=p/'external_dependencies.json';j=json.loads(file.read_text());j['average_ordered_neighbor_bound']=372;script='verify_all.py'
   elif name=='negative_count_coefficient':file=p/'radial/certificate.json';j=json.loads(file.read_text());j['lambda_count_exact']='-1';script='radial/verify.py'
   else:file=p/'radial/fit.json';j=json.loads(file.read_text());j['parameters'][0]+=1;script='verify_all.py'
   file.write_text(json.dumps(j));cmd=[sys.executable,str(p/script)]
  r=subprocess.run(cmd,capture_output=True,text=True,timeout=30);assert r.returncode!=0,(name,'corruption accepted');results[name]={'rejected':True,'returncode':r.returncode}
(D/'NEGATIVE_CONTROLS.json').write_text(json.dumps(results,indent=2)+'\n');print(results)
