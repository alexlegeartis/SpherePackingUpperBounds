from pathlib import Path
from fractions import Fraction as F
import json,math,sys
sys.set_int_max_str_digits(0);D=Path(__file__).parent;src=Path('research/lattice48_compact12_26_20261117');z=json.load(open(src/'harmonic_matrix_narrow_candidate.json'));co=list(map(F,next(r for r in json.load(open(src/'restrictions.json'))if r['n']==48)['coefficients'][::2]));q=sum(co);eps=F(1,10**8);guard=F(sys.argv[1])if len(sys.argv)>1 else F(1,10000)
Bs=[]
for B in z['B']:
 a=F(str(B[0][0]))+eps;b=(F(str(B[0][1]))+F(str(B[1][0])))/2;c=F(str(B[1][1]))+eps;assert a>0 and c>0 and a*c-b*b>0;Bs.append([[a,b],[b,c]])
Bs[0][0][0]-=guard;aa=Bs[0][0][0];bb=Bs[0][0][1];cc=Bs[0][1][1];det=aa*cc-bb*bb;assert aa>0 and det>0
zz=list(map(lambda x:F(str(x)),z['z']));zz[0]-=2*eps+F(1,10**9);ss=(cc*zz[0]**2-2*bb*zz[0]*zz[1]+aa*zz[1]**2)/det+F(1,10**12)
out=dict(status='PENDING UNIVERSAL GEOMETRY AND SCALAR REPLAY',co=list(map(str,co)),q=str(q),B=[[[str(x)for x in row]for row in B]for B in Bs],z=list(map(str,zz)),s=str(ss),eps=str(eps),constant_guard=str(guard),powers=[12,24],c=z['c'],degree=71,d=22)
(D/'frozen.json').write_text(json.dumps(out,indent=2));(D/'numerical_source.json').write_text(json.dumps(z,indent=2));print('PSD repaired; preview density',(sum(z['c'])+float(ss))*(22/4)**24/math.factorial(24),flush=True)
