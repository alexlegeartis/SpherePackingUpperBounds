import sys,json,time,math
sys.path.insert(0,'solver_packages');sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from flint import fmpq_poly,fmpq,arb_poly,arb,ctx
ctx.prec=384;p='continuum_three_point/';j=json.load(open(p+'adaptive_cyclic40_radial111.json'));P=fmpq_poly(list(map(fmpq,j['physical_polynomial'])));H=fmpq_poly(list(map(fmpq,j['fourier_polynomial'])));norm=H[0];P/=norm;H/=norm;B=fmpq(j['B']);M=fmpq(j['neighbor_cost'])/(B*norm)*fmpq(1000001,1000000);k=fmpq(F(str(j['kappa'])).numerator,F(str(j['kappa'])).denominator)
def atan(q):
 x=fmpq(q);lo=sum(((-1)**i*x**(2*i+1)/(2*i+1)for i in range(100)),fmpq(0));return lo,lo+x**201/201
al,au=atan('1/5');bl,bu=atan('1/239');pl=16*al-4*bu;pu=16*au-4*bl
assert 6*k*pl>40
for name,Q,cut in [('P',P,7*pl),('H',H,fmpq(0))]:
 roots=Q.numer().complex_roots();assert sum(m for z,m in roots)==Q.degree();assert all(not z.imag.contains(0)or z.real<arb(cut)for z,m in roots),name
assert P(7*pl)<0 and H[0]==1
AP=arb_poly([arb(c)for c in P.coeffs()]);stack=[(6*pl,7*pu,0)];passed=0;maxdepth=0;start=time.time();bad=[]
while stack:
 a,b,depth=stack.pop();c=(a+b)/2;rad=(b-a)/2;T=AP(arb_poly([arb(c),arb(rad)]));upper=T[0]+sum((abs(T[i])for i in range(1,len(T))),arb(0))
 if upper<0:passed+=1;maxdepth=max(maxdepth,depth);continue
 e=(((k-1)*arb(b)-6*k*arb(pl))).exp();cap=2-(-k*(arb(a)-6*arb(pu))).exp();v=upper*e-arb(M)*cap
 if v<0:passed+=1;maxdepth=max(maxdepth,depth);continue
 if depth>=32:
  cv=AP(arb(c))*((k-1)*arb(c)-6*k*arb(pl)).exp()-arb(M)*(2-(-k*(arb(c)-6*arb(pu))).exp())
  bad.append({'x':str(c),'value':str(cv)});print('FAIL',depth,str(cv),flush=True);break
 stack.extend([(a,c,depth+1),(c,b,depth+1)])
 if passed and passed%1000==0:print('cap',passed,len(stack),depth,time.time()-start,flush=True)
U=P[0]+B*M
out={'pass':not bad,'radial_file':'adaptive_cyclic40_radial111.json','kappa':str(k),'B_exact':str(B),'M_exact':str(M),'objective_exact':str(U),'objective_float':float(U),'pi_lower':str(pl),'pi_upper':str(pu),'global_roots_pass':True,'cap_rectangles':passed,'max_depth':maxdepth,'cap_failures':bad,'precision':ctx.prec,'elapsed':time.time()-start}
json.dump(out,open(p+'cyclic40_radial111_exact.json','w'),indent=2);print('RESULT',out['pass'],float(U),passed,maxdepth,time.time()-start,flush=True)
