from pathlib import Path
from fractions import Fraction as F
import random,itertools,math,subprocess,json,struct
def bits(x):return str(struct.unpack('>Q',struct.pack('>d',x))[0])
p=Path(__file__).parent;rng=random.Random(20261117)
def enc(q):
 v=float(q);l=math.nextafter(v,-math.inf);h=math.nextafter(v,math.inf);assert F.from_float(l)<=q<=F.from_float(h);return l,h
for trial in range(40):
 n=3;A=[[rng.randint(-2,2) for j in range(n)]for i in range(n)];G=[[sum(A[i][k]*A[j][k] for k in range(n))+int(i==j) for j in range(n)]for i in range(n)];r=[rng.randrange(2)for i in range(n)];cut=rng.randrange(30);mu=[[F(0)for j in range(n)]for i in range(n)];d=[]
 for i in range(n):
  d.append(F(G[i][i])-sum(mu[i][k]**2*d[k]for k in range(i)));mu[i][i]=1
  for j in range(i+1,n):mu[j][i]=(F(G[j][i])-sum(mu[j][k]*mu[i][k]*d[k]for k in range(i)))/d[i]
 expected=False
 for z in itertools.product(range(-5,6),repeat=n):
  v=[r[i]-2*z[i]for i in range(n)];q=sum(v[i]*G[i][j]*v[j]for i in range(n)for j in range(n))
  if q<=cut:expected=True;break
 lines=[str(n),str(cut),' '.join(map(str,r))]+[' '.join(map(str,row))for row in G]+[' '.join(map(bits,enc(v)))for v in d]+[' '.join(map(bits,enc(v)))for row in mu for v in row]
 fn=p/'test_input.txt';fn.write_text('\n'.join(lines)+'\n');out=subprocess.run([str(p/'verify_order2'),str(fn),'100000','fast','bits'],capture_output=True,text=True);a=json.loads(out.stdout);assert a['status'] in ['EXACT_WITNESS','EXHAUSTIVE_NO_VECTOR'];assert(a['status']=='EXACT_WITNESS')==expected
 if expected:
  v=a['residual_numerator'];assert all((v[i]-r[i])%2==0 for i in range(n));assert sum(v[i]*G[i][j]*v[j]for i in range(n)for j in range(n))==a['norm_numerator_over4']<=cut
print('PASS 40 exact brute-force small-instance comparisons')
