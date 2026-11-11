exec(open(__file__.replace('gaussian_screen','transfer_screen')).read().split('res=[]')[0])
import mpmath as mp
from scipy.optimize import minimize_scalar
mp.mp.dps=150
pm=[mp.mpf(x.numerator)/x.denominator for x in p];Mm=mp.mpf(M.numerator)/M.denominator
res=[]
for n in range(25,49):
 a=mp.mpf(n)/2;mom=F(1);b0=F(0)
 for j in range((len(pc)+1)//2):b0+=pc[2*j]*mom;mom*=F(2*j+1,n+2*j)
 Bn=mp.mpf((B/b0).numerator)/(B/b0).denominator
 rising=[mp.rf(a,k) for k in range(len(p))]
 def obj(t):
  w=1+mp.mpf(t);N=mp.fsum(x*y/w**k for k,(x,y) in enumerate(zip(pm,rising)))/w**a
  return (pm[0]+Mm*Bn*mp.exp(-6*mp.pi*t))/N*mp.pi**a*mp.mpf('1.5')**a/mp.gamma(a+1)
 fit=minimize_scalar(lambda t:float(obj(t)),bounds=(0,2),method='bounded',options={'xatol':1e-11})
 d=float(obj(fit.x));row={'n':n,'gaussian_parameter':fit.x,'density':d,'reference':bench[n],'ratio':d/bench[n]};res.append(row);print(row,flush=True)
(D/'gaussian_screen.json').write_text(json.dumps(res,indent=2))
