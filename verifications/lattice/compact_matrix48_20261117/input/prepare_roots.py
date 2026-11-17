import sys;sys.set_int_max_str_digits(0)
import json,math,pathlib,shutil
from fractions import Fraction as F
from flint import fmpq_poly as P,fmpq,ctx
n=int(sys.argv[1]);alpha=F(n,2)-1;ctx.prec=768;D=pathlib.Path('research/lattice48_compact_matrix_20261117/root_repair');D.mkdir(exist_ok=True);shutil.copy('research/lattice48_compact_matrix_20261117/frozen.json',D/'input_candidate.json');z=json.load(open(D/'input_candidate.json'));c=list(map(lambda x:F(str(x)),z['c']));x=P([0,1]);ls=[P([1]),P([fmpq(str(alpha+1)),-2])]
for k in range(2,len(c)):ls.append(((fmpq(str(2*k+alpha-1))-2*x)*ls[-1]-fmpq(str(k+alpha-1))*ls[-2])/k)
q=P([0]);norm=F(1)
for k,ck in enumerate(c):
 if k:norm*=F(k+alpha,k)
 q+=ls[k]*fmpq(str(ck*(-1)**k/norm))
json.dump({'Q':list(map(str,q))},open(D/'input_roots.json','w'))
s=pathlib.Path('research/lattice50_root_repair_20261117/repair.py').read_text().replace('Poly([25,-1])','Poly([Q("'+str(alpha+1)+'"),-1])').replace('2*k+23-x','2*k+Q("'+str(alpha-1)+'")-x').replace('(k+23)*Ls[-2]','(k+Q("'+str(alpha-1)+'"))*Ls[-2]').replace('dimension=50','dimension='+str(n));s=s.replace('len(neg)==1 and len(pos)%2==0 and 1+len(pos)+2*len(complexes)==old.degree()','len(neg)>=1 and len(pos)%2==0 and len(neg)+len(pos)+2*len(complexes)==old.degree()').replace("a=-rat(neg[0]);assert a>0;q=Poly([1,1/a]);factors=[dict(type='negative_linear',a=str(a))]","q=Poly([1]);factors=[]\nfor nr in neg:\n a=-rat(nr);assert a>0;q*=Poly([1,1/a]);factors.append(dict(type='negative_linear',a=str(a)))");(D/'repair.py').write_text(s)
