
if not __debug__: raise RuntimeError("Verification must not run under python -O.")
import json,sys,math
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_mat as Mat, fmpq_poly as Poly,arb,ctx
ctx.prec=256
from gaussian_guards import islands,outward,PI_LO,bernstein_upper,verify_guards
x=json.load(open(sys.argv[1]))[0];outdir=Path(sys.argv[2]);outdir.mkdir(exist_ok=True)
assert Q(str(F(x['exclusion_squared'])))==1
n=x['dimension'];degree=x['radial_degree'];assert degree==31 and len(x['coefficients'])==32 and n in [5,6,7], 'This exporter is restricted to the cited low-dimensional count tables';assert all(Q(str(F(str(v))))>=0 for v in x['count_weights']);R=Q(str(F(x['truncation_squared'])));r=F(str(math.sqrt(float(R))));assert Q(str(r*r))==R
M=[[[str(F(str(v)))for v in row]for row in B]for B in x['matrices']]
# Check every exact block.
checks=[]
for ell,B in enumerate(M):
 A=Mat([[Q(c)for c in row]for row in B]);assert A==A.transpose();N=A.nrows();assert N==x['radial_nodes']+int(ell==0);L=Mat(N,N);d=[]
 for j in range(N):
  L[j,j]=1;t=A[j,j]-sum((L[j,k]**2*d[k]for k in range(j)),Q(0));d.append(t);assert t>0,(ell,j)
  for i in range(j+1,N):L[i,j]=(A[i,j]-sum((L[i,k]*L[j,k]*d[k]for k in range(j)),Q(0)))/t
 checks.append({'block':ell,'positive_definite':True,'min_pivot':float(min(d))})
X=Poly([0,1]);a=Q(n,2);ls=[Poly([1]),a-2*X]
for k in range(1,degree):ls.append(((2*k+a-2*X)*ls[-1]-(k+a-1)*ls[-2])/Q(k+1))
c=list(map(lambda v:Q(str(F(str(v)))),x['coefficients']));P=sum((c[k]*ls[k]/ls[k][0]for k in range(degree+1)),Poly([]));H=sum(((-1)**k*c[k]*ls[k]/ls[k][0]for k in range(degree+1)),Poly([]))
guardc=Q('314/100');phys=islands(P,R*PI_LO,-1);y0=(2*n-guardc)/16;freq=islands(H,y0,1);eta=zeta=Q(0)
for lo,hi in phys:eta=max(eta,outward(arb(bernstein_upper(P,lo,hi))*(-arb(3)*arb(lo)/4).exp()/arb(lo-guardc)))
for lo,hi in freq:zeta=max(zeta,outward(arb(2**n)*arb(bernstein_upper(-H,lo,hi))*(-arb(3)*arb(lo)/4).exp()/arb(guardc-Q(n,8)+lo/16)))
D=Q(2**n)*(2*n-guardc);E=(guardc-Q(n,8))/Q(2**n);norm=H[0]-eta*D+zeta*E
guards={'gaussian_guard_family':'reciprocal_gaussians','dimension':n,'reciprocal_gaussian_guards':{'c':str(guardc),'eta':str(eta),'zeta':str(zeta),'physical_intervals':[[str(a),str(b)]for a,b in phys],'fourier_intervals':[[str(a),str(b)]for a,b in freq],'corrected_normalization':str(norm)}}
verify_guards(guards,P,H,Q(1),R,PI_LO)
count_radii=[str(F(str(math.sqrt(z))).limit_denominator(10**6))for z in x['count_squared_radii']]
tables={5:([F('1.003012'),F('1.008096'),F('1.013258'),F('1.017445')],[45,46,47,48]),6:([F('1.002005'),F('1.004023'),F('1.006053'),F('1.008096')],[79,80,81,82]),7:([F('1.001000'),F('1.002005'),F('1.003012'),F('1.004023'),F('1.005036')],[136,137,138,139,140])}
if n in tables and x['count_weights']:
    assert list(map(F,count_radii))==tables[n][0] and x['count_bounds']==tables[n][1]
source={'kernel_family':'chebyshev_radius_count_cyclic','dimension':n,'R':str(r),'J':1,'p':x['radial_nodes']-1,'K':x['angular_degree'],'x':x['coefficients'],'count_radii':count_radii,'count_weights':x['count_weights'],'count_bounds':x['count_bounds'],'scope':'UNVERIFIED Chebyshev cyclic/count candidate; exact prerequisites only.'}
sourcepath=outdir/'witness.json';sourcepath.write_text(json.dumps(source,indent=2))
(outdir/'witness_psd_exact.json').write_text(json.dumps({'source':'witness.json','matrices':M,'checks':checks},indent=2))
B=Mat([[Q(c)for c in row[1:]]for row in M[0][1:]]);b=Mat([[Q(row[0])]for row in M[0][1:]]);alpha=Q(M[0][0][0])-(b.transpose()*B.inv()*b)[0,0];assert alpha>=0
count_cost=sum((Q(str(F(str(l))))*Q(str(F(str(b))))for l,b in zip(x['count_weights'],x['count_bounds'])),Q(0))
raw=(P[0]+Q(M[0][0][0])+count_cost+guardc*(eta+zeta))/norm;saturated=raw-alpha/norm;volume=arb.pi()**(arb(n)/2)/(arb(n)/2+1).gamma()/arb(2)**n
out={'scope':'Unverified until complete pair/far/cyclic checks pass','dimension':n,'guards':guards,'raw_density':str(volume*arb(raw)),'saturated_density':str(volume*arb(saturated)),'origin_slack':str(alpha),'raw_objective':str(raw),'saturated_objective':str(saturated)}
(outdir/'prerequisites.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
