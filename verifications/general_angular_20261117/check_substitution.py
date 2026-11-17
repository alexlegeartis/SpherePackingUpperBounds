from fractions import Fraction as F
from pathlib import Path
from math import factorial
import json
P=Path(__file__).resolve().parent
S=P/'source'
def atan(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N)),F(0))
piU=16*atan(5,201)-4*atan(239,200)
for n in [45,48]:
 new=json.loads((P/f'dim{n}/certificate.json').read_text());old=json.loads((S/f'certificate{n}.json').read_text())
 for k in ['physical_polynomial','fourier_polynomial','s_exact','a_exact','b_exact','q_exact','cut_X_exact']:assert new[k]==old[k],k
 C,B=F(new['C_exact']),F(new['B_exact']);assert C/B==F(old['C_exact'])/F(old['B_exact']) and C<F(old['C_exact'])
 p0=F(new['physical_polynomial'][0]);h0=F(new['fourier_polynomial'][0]);s=F(new['s_exact']);z=(p0+C)/h0;bound=F(new['strict_density_upper_bound'])
 assert z==F(new['objective_exact'])>0 and 0<bound<F(old['strict_density_upper_bound'])
 if n%2==0:assert z*(piU*s/4)**(n//2)/factorial(n//2)<bound
 else:assert z*z*4*piU**(n-1)*s**n*F(factorial((n+1)//2),factorial(n+1))**2<bound*bound
 print(n,'SOURCE SUBSTITUTION AND INDEPENDENT201/200-TERM MACHIN ROUNDING PASS',str(bound))
