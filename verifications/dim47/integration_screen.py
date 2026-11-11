exec(open(__file__.replace('integration_screen','transfer_screen')).read().split('res=[]')[0])
import mpmath as mp
from scipy.optimize import minimize_scalar
mp.mp.dps=130
pm=[mp.mpf(x.numerator)/x.denominator for x in p];hm=[mp.mpf(F(x).numerator)/F(x).denominator for x in r['fourier_polynomial']];Mm=mp.mpf(M.numerator)/M.denominator;kappa=mp.mpf(17)/8
res=[]
for n in range(25,48):
 a=mp.mpf(n)/2;b=mp.mpf(48-n)/2;mom=F(1);b0=F(0)
 for j in range((len(pc)+1)//2):b0+=pc[2*j]*mom;mom*=F(2*j+1,n+2*j)
 Bn=mp.mpf((B/b0).numerator)/(B/b0).denominator
 rising=[mp.rf(b,k) for k in range(len(p))]
 def obj(t):
  t=mp.mpf(t);w=1+t
  f0=mp.fsum(x*y/w**k for k,(x,y) in enumerate(zip(pm,rising)))/w**b
  N=mp.fsum(x*y*(t/w)**k for k,(x,y) in enumerate(zip(hm,rising)))/w**b
  cost=Mm*Bn*(2*kappa+t)**b/(kappa+t)**(2*b)
  return (f0+cost)/N*mp.pi**a*mp.mpf('1.5')**a/mp.gamma(a+1)
 fit=minimize_scalar(lambda t:float(obj(t)),bounds=(0,100),method='bounded',options={'xatol':1e-9})
 vals=[(float(obj(0)),0),(float(obj(fit.x)),float(fit.x))];d,t=min(vals)
 row={'n':n,'gaussian_parameter':t,'density':d,'reference':bench[n],'ratio':d/bench[n]};res.append(row);print(row,flush=True)
(D/'integration_screen.json').write_text(json.dumps(res,indent=2))
