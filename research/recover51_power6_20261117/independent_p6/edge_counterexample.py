from pathlib import Path
import json,hashlib
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=512
D=Path(__file__).parent;K=D/'rejected_kernel.json';v=json.load(open(K));assert v['powers']==[6,10,14] and v['R']=='19/13' and v['dimension']==51
ar=lambda x:arb(F(x).numerator)/F(x).denominator
co=list(map(ar,v['co']));B=[[[ar(x)for x in row]for row in b]for b in v['B']]
a=arb(1);b=arb(1);c=1+ar('6/13')*ar('5/8');beta=(1,1,ar('3/8'))
def angular(t):
 G=[arb(1),t]
 for k in range(2,30):G.append(((2*k+51-4)*t*G[-1]-(k-1)*G[-2])/(k+51-3))
 return G
def term(t,u,w):
 G=angular(t);U=[u**p for p in [6,10,14]];W=[w**p for p in [6,10,14]]
 return sum(co[k]*G[k]*U[i]*B[k][i][j]*W[j]for k in range(30)for i in range(3)for j in range(3))
value=term((a+b-c)/(2*(a*b).sqrt()),beta[0],beta[1])+term((a+c-b)/(2*(a*c).sqrt()),beta[0],beta[2])+term((b+c-a)/(2*(b*c).sqrt()),beta[1],beta[2])
assert value>0
out={'status':'REJECTED_POSITIVE_CYCLIC_VALUE','cube_point':['0','0','5/8'],'squared_sides':['1','1',str(F(1)+F(6,13)*F(5,8))],'raw_cyclic_value':str(value),'precision':512,'frozen_sha256':hashlib.sha256(K.read_bytes()).hexdigest(),'method':'independent direct nine-entry feature contraction, no Jet or gate code'}
(D/'edge_counterexample.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
