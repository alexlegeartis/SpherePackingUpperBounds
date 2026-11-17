from pathlib import Path
from fractions import Fraction as F
import random,math,struct,itertools,subprocess,json
P=Path(__file__).resolve().parent;random.seed(20261117);out=[]
def bits(x):return struct.unpack('Q',struct.pack('d',x))[0]
def enclosure(x):
 y=float(x);lo=math.nextafter(y,-math.inf);hi=math.nextafter(y,math.inf);assert F(lo)<=x<=F(hi);return f'{bits(lo)} {bits(hi)}'
for test in range(40):
 n=1+test%4;scale=[1,8,40,100][test//4%4];A=[[random.randrange(-2,3)for j in range(n)]for i in range(n)];G=[[scale*(sum(A[k][i]*A[k][j]for k in range(n))+(i==j))for j in range(n)]for i in range(n)];r=[random.randrange(2)for i in range(n)]
 L=[[F(i==j)for j in range(n)]for i in range(n)];d=[]
 for i in range(n):
  d.append(F(G[i][i])-sum(L[i][k]**2*d[k]for k in range(i)))
  for j in range(i+1,n):L[j][i]=(F(G[j][i])-sum(L[j][k]*L[i][k]*d[k]for k in range(i)))/d[i]
 expected=any(sum((r[i]-2*z[i])*G[i][j]*(r[j]-2*z[j])for i in range(n)for j in range(n))<=24 for z in itertools.product(range(-4,5),repeat=n))
 p=P/'control_input.txt';p.write_text(f'{n} 24\n'+' '.join(map(str,r))+'\n'+'\n'.join(' '.join(map(str,row))for row in G)+'\n'+'\n'.join(enclosure(x)for x in d)+'\n'+'\n'.join(enclosure(x)for row in L for x in row)+'\n')
 q=subprocess.run([str(P/'controls_search'),str(p)],text=True,capture_output=True);assert q.returncode in [0,10],q.stdout;status=json.loads(q.stdout.splitlines()[-1])['status'];assert (status=='VECTOR')==expected,(test,q.stdout,expected);out.append(dict(test=test,dimension=n,expected_vector=expected,status=status))
(P/'controls_result.json').write_text(json.dumps(out,indent=2));print('PASS40 exact brute-force controls',sum(r['expected_vector']for r in out),'nonempty')
