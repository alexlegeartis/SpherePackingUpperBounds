from pathlib import Path
import json,math
from fractions import Fraction as F
p=Path(__file__).resolve().parent.parent
old=json.loads((p/'certificate48.json').read_text());c=json.loads((p/'combined48_trial/certificate.json').read_text());P=list(map(F,c['p']));deg=len(P)-1;rise=[math.prod(range(24,24+j))for j in range(deg+1)];D=[F(0)]*(deg+1)
for j in range(deg,-1,-1):
 D[j]=((-1)**j*P[j]*math.factorial(j)*rise[j]-sum(D[k]*F(math.factorial(k),math.factorial(k-j))for k in range(j+1,deg+1)))/math.factorial(j)
Q=[(-1)**j*sum(D[k]*(-1)**k*F(math.factorial(k),math.factorial(k-j))for k in range(j,deg+1))/(math.factorial(j)*rise[j])for j in range(deg+1)]
oldq=list(map(F,old['q']));tau=F(c['tau']);v=F(c['v']);w=F(c['w'])
expected=[sum(m*(oldq[j-k]if 0<=j-k<len(oldq)else 0)for k,m in enumerate((F(c['fourier_multiplier'][0]),v/tau,w/tau**2)))for j in range(deg+1)]
m=list(map(F,c['fourier_additive_coefficients']));expected=[a+F(c['qzero'])*(m[j]/tau**j if j<len(m)else 0)for j,a in enumerate(expected)]
assert Q==expected and Q[0]==F(c['qzero'])
print('PASS independent exact Laguerre transform, degree',deg,'and unchanged Fourier zero normalization')
