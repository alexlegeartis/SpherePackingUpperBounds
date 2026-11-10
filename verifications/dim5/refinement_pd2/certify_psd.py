import sys,json,math
from fractions import Fraction as F
from decimal import Decimal,localcontext
r=json.load(open(sys.argv[1]));off=32;out=[];matrices=[]
with localcontext()as ctx:
 ctx.prec=90;rt=Decimal(2).sqrt();den=10**30
 for k in range(r['K']+1):
  n=r['J']*(r['p']+1)+(k==0);A=[[F(0)for _ in range(n)]for _ in range(n)]
  for j in range(n):
   for i in range(j+1):
    z=Decimal(str(r['x'][off]))/(1 if i==j else rt);off+=1;A[i][j]=A[j][i]=F(int(z*den),den)
  matrices.append([[str(v)for v in row]for row in A]);d=[];L=[[F(int(i==j))for j in range(n)]for i in range(n)];ok=True
  for j in range(n):
   t=A[j][j]-sum(L[j][h]**2*d[h]for h in range(j));d.append(t)
   if t<=0:ok=False;break
   for i in range(j+1,n):L[i][j]=(A[i][j]-sum(L[i][h]*L[j][h]*d[h]for h in range(j)))/t
  out.append({'block':k,'positive_definite':ok,'min_pivot':float(min(d))})
print(json.dumps(out),flush=True);json.dump({'source':sys.argv[1],'checks':out,'matrices':matrices,'scope':'Exact PSD of rationalized matrices only; kernel signs not certified'},open(sys.argv[1].replace('.json','_psd_exact.json'),'w'),indent=2)
