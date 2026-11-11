from pathlib import Path
import json,hashlib
from quantify_schur import audit
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent;root=p.parent.parent
for name,digest in json.loads((p/'source_inputs.json').read_text()).items():
 assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
registry={r['dimension']:r for r in json.loads((root/'RESULTS.json').read_text())}
results=[]
for n,bound in {18:'0.013224473',19:'0.009521000',20:'0.006936042'}.items():
 source=root/f'verifications/dim{n}/refinement_cyclic/certificate.json'
 result=audit(source,bound);assert result['dimension']==n
 assert registry[n]['strict_density_upper_bound']==bound and registry[n]['package']=='verifications/cyclic_origin_slack'
 result['source']=str(source.relative_to(root));result['registration_status']='Exact saturation verified; run verify.py to replay all source domains'
 results.append(result);print(n,result['density_interval'],bound,flush=True)
(p/'saturation_results.json').write_text(json.dumps(results,indent=2)+'\n')
import subprocess,sys
subprocess.run([sys.executable,str(p/'verify_rational.py')],check=True)
print('PASS exact origin saturation, independent rational rounding and all three strict outward thresholds')
