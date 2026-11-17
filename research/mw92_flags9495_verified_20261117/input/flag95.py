from pathlib import Path
from fractions import Fraction as F
import json,math,mpmath as mp
p=Path(__file__).parent;r=json.load(open(p/'results.json'));h=[16-F(r['depths'][str(j)])for j in[1,3,9]];a,b,c=h;det=a*b*c-a*b*b/4-a*a*c/4;center2=F(2**184,47**14)*64/det;mp.mp.dps=65;fm=lambda x:mp.mpf(x.numerator)/x.denominator;occ=mp.sqrt(fm(center2))*mp.pi**mp.mpf('47.5')/mp.gamma(mp.mpf('48.5'));print('density95',occ);pi=F(r['pi_upper']);atan=lambda q,N:sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0));lo=16*atan(5,100)-4*atan(239,101);df=math.prod(range(1,96,2));floor=F('2.06770140277194e-21');sq=center2*F(2**96,df**2)*lo**94
floor=F(str(mp.floor(occ*mp.mpf(10)**35)))/10**35
assert floor**2<sq
print('floor',floor,float(floor));old=F('2.043145398811e-21');print('gain%',100*(occ/fm(old)-1));json.dump({'gram':[[str(a),str(-a/2),'0'],[str(-a/2),str(b),str(-b/2)],['0',str(-b/2),str(c)]],'determinant':str(det),'center_density_squared':str(center2),'occupied95_floor':str(floor),'occupied95_decimal':str(occ)},open(p/'flag95.json','w'),indent=2)
