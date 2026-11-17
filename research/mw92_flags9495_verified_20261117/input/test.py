from fractions import Fraction as F
import math,json,mpmath as mp
mp.mp.dps=60;n=92;a=F(n,2);mu=F(2,n+1);m2=F(24*n+4*n*(n-1),math.prod(range(n,n+4)));m3=F(720*n+144*n*(n-1)+8*n*(n-1)*(n-2),math.prod(range(n,n+6)));v=m2-mu*mu;c=m3-3*mu*m2+2*mu**3
C=(n+1)*F(n,n+3)**(n//2)*(1+a*(a+1)*v/(2*(1+mu)**2)-a*(a+1)*(a+2)*c/(6*(1+mu)**3))
atan=lambda q,N:sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
pi_upper=16*atan(5,101)-4*atan(239,100);V=pi_upper**46/math.factorial(46);to_mp=lambda x:mp.mpf(x.numerator)/x.denominator
D={}
for j in[1,2,3,9,27]:
 value=(to_mp(C)*47**7/(j*to_mp(V)))**(mp.mpf(1)/46);d=F(int(value*10**9),10**9)
 while j*V*d**46>=C*47**7:d-=F(1,10**9)
 D[j]=d
h1=16-D[1];h2=16-D[3];det=h1*h2-h1*h1/4;old=F(3,4)*h2*h2
assert h1>4 and h2>h1 and h2<2*h1
center2=F(2**184,47**14)*16/det
occupied=mp.sqrt(to_mp(center2))*mp.pi**47/mp.factorial(47)
print('h1,h2',float(h1),float(h2),'gain%',100*(mp.sqrt(to_mp(old/det))-1),'occupied94',occupied)
json.dump({'covering_constant':str(C),'pi_upper':str(pi_upper),'depths':{str(k):str(v)for k,v in D.items()},'gram':[[str(h1),str(-h1/2)],[str(-h1/2),str(h2)]],'determinant':str(det),'center_density_squared':str(center2),'occupied94':str(occupied),'gain_over_uniform_parity_percent':str(100*(mp.sqrt(to_mp(old/det))-1))},open('research/mw92_parallelogram_layers_20261117/results.json','w'),indent=2)
