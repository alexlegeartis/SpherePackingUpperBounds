"""Prepare a NEW exact producer check; never inherit historical pass receipts."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,shutil,itertools
import sympy as s
D=Path(__file__).parent;candidate=Path(sys.argv[1]);v=json.load(open(candidate));assert v['status']=='NUMERICAL_ONLY_NOT_A_BOUND' and v['minimum_eigenvalue']>0
p=v['powers'][0];assert v['powers']==[p,p+4,p+8]
out=D/(sys.argv[2] if len(sys.argv)>2 else f'producer_p{p}');out.mkdir(exist_ok=True)
base=D/'base';ker=json.load(open(base/'geometry/dim51/frozen.json'));ker['powers']=v['powers'];rat=lambda x:str(F(format(float(x),'.16g')))
ker['B']=[[[rat((b[i][j]+b[j][i])/2)for j in range(3)]for i in range(3)]for b in v['B']];ker['z']=list(map(rat,v['z']))
for b in ker['B']:
 m=s.Matrix(b)
 for k in range(1,4):
  for ids in itertools.combinations(range(3),k):assert m.extract(ids,ids).det()>0
z=s.Matrix(ker['z']);origin=(z.T*s.Matrix(ker['B'][0]).inv()*z)[0];ss=F(str(origin))+F(1,10**9);ker['s']=str(ss)
G=out/'geometry';S=out/'scalar';(G/'dim51').mkdir(parents=True,exist_ok=True);S.mkdir(exist_ok=True)
for name in ['jet.py','jet_fast.py','verify_caps.py','verify_edges.py','gate.py']:
 src=(base/'geometry'/name).read_text();src=src.replace('[8,12,16]',str(v['powers'])).replace('(1-z)**8','(1-z)**'+str(p)).replace('((1-y)*(1-z))**8','((1-y)*(1-z))**'+str(p));(G/name).write_text(src)
(G/'dim51/frozen.json').write_text(json.dumps(ker,indent=2));(S/'kernel.json').write_text(json.dumps(ker,indent=2));shutil.copy(base/'scalar/repaired.json',S/'repaired.json')
for name in ['certify_annuli.py','verify_fourier_tail.py','build_caps.py','verify_caps.py']:
 src=(base/'scalar'/name).read_text().replace('[8,12,16]',str(v['powers']))
 if name=='certify_annuli.py':src=src.replace('**4\nkp=', '**1\nkp=').replace('h**(i+2)',f'h**({p}+4*i)').replace('h**(i+j+4)',f'h**({2*p}+4*i+4*j)')
 (S/name).write_text(src)
(out/'STATUS.md').write_text('New reconstructed producer candidate. Exact block PSD checked. All geometry, annular, tail, and independent audits pending. NOT A BOUND.\n')
(out/'source.json').write_text(json.dumps({'candidate':candidate.name,'source_commit':'34f59648d034e92a52768c95af96221a8f589ac9','s':ker['s']},indent=2))
print(out)
