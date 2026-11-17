from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import re,html,json,hashlib
from flint import fmpz_mat as Z
P=Path(__file__).resolve().parent;source=P/'L80_source.html';s=source.read_text();body=s.split('<a NAME="GRAM"',1)[1].split('</a><br>',1)[1].split('<a NAME=',1)[0];v=list(map(int,re.findall(r'-?\d+',html.unescape(re.sub('<[^>]*>',' ',body)))));assert v[:2]==[80,0]and len(v)==2+80*81//2
G=[[0]*80 for _ in range(80)];pos=2
for i in range(80):
 for j in range(i+1):G[i][j]=G[j][i]=v[pos];pos+=1
M=Z(G);assert M.det()==1 and all(Z([r[:k]for r in G[:k]]).det()>0 for k in range(1,81))
assert [r[:2]for r in G[:2]]==[[8,-4],[-4,8]]
# First two coordinate basis vectors span a primitive2A2 normal sublattice.
K=Z([r[:2]for r in G[:2]]);assert K.det()==48
old=F(2**78,48);assert old==F(2**74,3)
new=F(json.loads((P/'result.json').read_text())['center_density_squared'])
with localcontext()as c:
 c.prec=70;gain=((D(new.numerator*old.denominator)/D(new.denominator*old.numerator)).sqrt()-1)*100
out={'status':'PASS exact pinned external-derived example; not strongest-prior assertion','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'catalogue_url':'https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html','published_parent_minimum':8,'ambient_det':1,'primitive_normal':'first two catalogue basis vectors, Gram[[8,-4],[-4,8]]','normal_det':48,'external_center_density_squared':str(old),'increase_percent_approx':str(gain)};(P/'external_result.json').write_text(json.dumps(out,indent=2));(P/'L80_source.html').write_bytes(source.read_bytes());print(json.dumps(out,indent=2))
