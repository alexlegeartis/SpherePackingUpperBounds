import sys,json,math
sys.path.insert(0,'solver_packages')
from flint import fmpq,fmpq_poly
from fractions import Fraction as F
r=json.load(open(sys.argv[1]));D=31;R=F(str(r['R']));X=fmpq_poly([0,1]);ls=[fmpq_poly([1]),2-2*X]
for k in range(1,D):ls.append(((2*k+2-2*X)*ls[-1]-(k+1)*ls[-2])/fmpq(k+1))
P=fmpq_poly([]);H=fmpq_poly([])
for k in range(D+1):
 c=fmpq(str(F(str(r['x'][k]))))/(k+1);P+=c*ls[k];H+=(-1)**k*c*ls[k]
def sturm(p):
 seq=[p,p.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);num=rem.numer();cont=abs(num.content());seq.append(fmpq_poly(num/cont))
 return seq
def var(v):
 ss=[1 if a>0 else -1 for a in v if a];return sum(a!=b for a,b in zip(ss,ss[1:]))
def roots(p):
 seq=sturm(p);return var([p[0]for p in seq])-var([p[p.degree()]for p in seq])
pi_lo=F(314159265358979323846264338327950288419716939937510,10**50);T=fmpq(str(R*R*pi_lo));out={'fourier_positive_roots':roots(H),'fourier_origin_positive':bool(H[0]>0),'tail_roots_after_conservative_threshold':roots(-P(X+T)),'tail_start_negative':bool(P(T)<0),'scope':'Exact rationalized coefficient polynomials; root counts must be zero plus signs correct to certify'}
print(json.dumps(out),flush=True);json.dump(out,open(sys.argv[1].replace('.json','_radial_exact.json'),'w'),indent=2)
