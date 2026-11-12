from math import comb,ceil
from fractions import Fraction as F
import scipy.optimize as so, numpy as np,json
from pathlib import Path
P=Path(__file__).resolve().parent
base=json.load(open(P.parent/'density_lower_baseline/records_1_48.json'))['general'];base={r['dimension']:r for r in base}
def K(n,k,x):return sum((-1)**j*3**(k-j)*comb(x,j)*comb(n-x,k-j)for j in range(max(0,k-(n-x)),min(k,x)+1))
rows=[]
for n in range(13,25):
 # dual polynomial 1+sum c_k K_k, f(i)<=0 for i>=8; nonnegative c_k
 # objective f(0) upper, normalize coefficient f0=1
 A=np.array([[K(n,k,i)for k in range(1,n+1)]for i in range(8,n+1)],dtype=float);c=np.array([K(n,k,0)for k in range(1,n+1)],dtype=float)
 # normalized basis K_k/K_k(0), objective sum variables, nonneg
 R=A/c
 r=so.linprog(np.ones(n),A_ub=R,b_ub=-np.ones(n-7),bounds=(0,None),method='highs')
 if not r.success:print(r.message);continue
 # rational coefficient reconstruction and rigorous repair by scaling if all negatives
 f=[F(float(max(0,x))).limit_denominator(10**9)/int(c[k])for k,x in enumerate(r.x)]
 vals=[sum(f[k-1]*K(n,k,i)for k in range(1,n+1))for i in range(8,n+1)]
 worst=max(vals);assert worst<0
 factor=max(F(1),-1/worst);f=[x*factor for x in f];u=1+sum(f[k-1]*K(n,k,0)for k in range(1,n+1));assert all(1+sum(f[k-1]*K(n,k,i)for k in range(1,n+1))<=0 for i in range(8,n+1))
 b2=F(base[2*n]['center_density_squared_exact']);target2=b2*3**(n+2);floor=int(np.sqrt(float(target2)))
 while F((floor+1)**2)<=target2:floor+=1
 while F(floor**2)>target2:floor-=1
 target=floor+1
 print(2*n,target,float(u), 'exclude' if u<target else 'open',flush=True)
 rows.append({'dimension':2*n,'length':n,'needed':target,'code_upper_rational':str(u),'coefficients':[str(x)for x in f],'excluded':u<target})
json.dump(rows,open(P/'code_sweep.json','w'),indent=2)
