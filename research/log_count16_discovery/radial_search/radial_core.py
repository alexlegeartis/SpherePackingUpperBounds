"""One fixed-profile degree-79 enrichment and matched degree-63 control.

Discovery and exact rejection diagnostics only. No angular refit and no
density certificate. Set BLAS/OMP threads to one when reproducing.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys, json, math, time, hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import eval_genlaguerre, gammaln
from flint import fmpq as Q, fmpq_poly as Poly, arb, ctx

sys.set_int_max_str_digits(0)
ctx.prec=320
D=Path(__file__).resolve().parent
I=D/'inputs'
COMPARATOR=.02492121
VK=np.array([1.02,1.05,1.1,1.2,1.4,2.])
LK=np.array([1.2,1.4,1.5,1.6,1.65,1.75])
PARS=np.array(json.loads((I/'fit.json').read_text())['parameters'])
B=float(Fraction(json.loads((I/'angular_dim16.json').read_text())['B_exact']))

def write(name,obj):
    (D/name).write_text(json.dumps(obj,indent=2)+'\n')

def rat(a):
    z=Fraction(float(a))
    return Q(z.numerator,z.denominator)

def basis(v):
    v=np.asarray(v)
    return np.concatenate([(v-1)[...,None],np.maximum(v[...,None]-VK,0),
         np.log(v)[...,None],np.maximum(np.log(v[...,None]/LK),0)],axis=-1)

def lag(x,d):
    norm=np.exp(gammaln(np.arange(d+1)+8)-gammaln(8)-gammaln(np.arange(d+1)+1))
    return np.array([eval_genlaguerre(k,7,2*np.asarray(x))/norm[k] for k in range(d+1)]).T

def polys(radial):
    d=len(radial)-1; X=Poly([0,1]); E=[Poly([1]),8-2*X]
    for k in range(1,d):
        E.append(((2*k+8-2*X)*E[-1]-(k+7)*E[-2])/(k+1))
    E=[p/p(0) for p in E]
    P=sum((rat(c)*e for c,e in zip(radial,E)),Poly())
    H=sum(((-1)**k*rat(c)*e for k,(c,e) in enumerate(zip(radial,E))),Poly())
    return P,H

def roots(poly,lo):
    all_roots=poly.numer().complex_roots()
    assert sum(m for z,m in all_roots)==poly.degree()
    result=[]
    for z,m in all_roots:
        if z.imag.contains(0) and not z.real<arb(lo):
            # Real roots from the isolator must be real, not unresolved boxes.
            assert z.imag.is_exact() and z.imag==0
            result.append(dict(x=float(z.real.mid()),multiplicity=m,
                real=str(z.real),imag=str(z.imag)))
    return sorted(result,key=lambda r:r['x'])

def audit(case,name):
    start=time.monotonic(); P,H=polys(case['radial'])
    write(name+'_exact.json',dict(physical=list(map(str,P.coeffs())),
        fourier=list(map(str,H.coeffs())),C=str(rat(case['C'])),
        normalization=str(H(0)),coefficient_rule='exact binary64 rational',certified_bound=False))
    out=dict(certified_bound=False,pair_proved=False,compact_cyclic_proved=False,gates={})
    for gate,p,lo,sign in [('fourier',H,0,1),('physical',P,40,-1)]:
        rr=roots(p,lo); witness=[]
        xs=[float(lo)]+[r['x'] for r in rr]
        for l,r in zip(xs,xs[1:]+[max(xs+[2500])+1]):
            x=rat((l+r)/2); val=p(x)
            if sign*val<0:
                witness.append(dict(x_exact=str(x),x=float(x),value_exact=str(val),
                     value_interval=str(arb(val)),function_value=str(arb(val)*(-arb(x)).exp())))
        passed=(not rr and sign*p(Q(lo))>0)
        out['gates'][gate]=dict(pass_conservative_root_gate=passed,
            roots=rr,negative_sign_interval_count=len(witness),counterexamples=witness)
        print('AUDIT',name,gate,'roots',len(rr),'bad intervals',len(witness),flush=True)
    out['elapsed_seconds']=time.monotonic()-start
    write(name+'_audit.json',out)
    return out,P,H

def shared_rows(source,P,H):
    original=np.load(I/'dim16_radial_rows.npz')
    extras={}; detail={}
    for gate,p,lo in [('fourier',H,0),('sign',P,40)]:
        critical=[]
        for kind,dp in [('polynomial_derivative',p.derivative()),
                        ('gaussian_derivative',p.derivative()-p)]:
            r=roots(dp,lo);detail[gate+'_'+kind]=r
            critical.extend(x['x'] for x in r)
        extras[gate]=np.array(critical)
    # Pair residual extrema on the finite annulus, numerically located only.
    c=np.array(source['radial']); C=source['C']
    def pairval(x):
        xx=np.atleast_1d(x);h=np.exp(-basis(xx/(2*np.pi))@PARS)
        return np.exp(-xx)*(lag(xx,63)@c)-C*(2*h-h*h)/B
    xx=np.unique(np.r_[np.linspace(2*np.pi,40,16001),2*np.pi*VK,2*np.pi*LK])
    vv=pairval(xx);inds=np.where((vv[1:-1]>=vv[:-2])&(vv[1:-1]>=vv[2:]))[0]+1
    cp=[]
    for i in inds:
        fit=minimize_scalar(lambda x:-float(pairval(x)[0]),bounds=(xx[i-1],xx[i+1]),
            method='bounded',options={'xatol':1e-13})
        cp.append(dict(x=float(fit.x),residual=float(-fit.fun)))
    extras['pair']=np.r_[[v['x'] for v in cp],2*np.pi*VK,2*np.pi*LK]
    detail['pair_local_maxima']=cp
    rows={k:np.unique(np.r_[original[k],extras[k]]) for k in ['pair','fourier','sign']}
    np.savez_compressed(D/'shared_rows.npz',**rows)
    write('source_critical_rows.json',dict(details=detail,counts={k:len(v) for k,v in rows.items()},
        original_counts={k:len(original[k]) for k in rows},
        description='Source critical roots and numerical pair maxima appended identically to both degrees.'))
    return rows

def solve(rows,d,s,name,source=None):
    start=time.monotonic();x,y,z=rows['pair'],rows['fourier'],rows['sign']
    h=np.exp(-basis(x/(np.pi*s))@PARS);g=np.exp(x)*(2*h-h*h)/B
    sg=(-1.)**np.arange(80)
    AA=np.vstack([np.c_[lag(x,79),-g],np.c_[lag(z,79),np.zeros(len(z))],
        np.c_[-lag(y,79)*sg,np.zeros(len(y))]])
    A=np.c_[AA[:,:d+1],AA[:,-1]]
    scale=np.maximum(np.abs(A).max(axis=1),1e-200)
    A/=scale[:,None]
    eq=np.r_[(-1.)**np.arange(d+1),0.]
    nested={}
    if source is not None:
        v=np.r_[source['radial'],np.zeros(d+1-len(source['radial'])),source['C']]
        nested=dict(padded_source_coefficients=v.tolist(),source_normalization=float(eq@v),
            source_objective=float(v.sum()),max_scaled_residual=float(np.max(A@v)),
            max_constraint_violation=float(np.max(A@v+1e-9)))
    np.savez_compressed(D/(name+'_lp_rows.npz'),pair=x,fourier=y,sign=z,
         scales=scale,A=A,b=np.full(len(A),-1e-9),equality=eq,objective=np.ones(d+2))
    fit=linprog(np.ones(d+2),A_ub=A,b_ub=np.full(len(A),-1e-9),
        A_eq=eq[None,:],b_eq=[1.],bounds=[(None,None)]*(d+1)+[(0,None)],method='highs',
        options={'time_limit':120,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    result=dict(certified_bound=False,success=bool(fit.success),message=fit.message,
        dimension=16,degree=d,s=s,cut_X=40.,row_count=len(A),source_nesting=nested,
        scaling='Per-degree maximum absolute row entry, exactly the source scaling rule; identical sample locations.',
        normalization='sum((-1)^k c_k)=1; normalized Laguerre L_k^7(2X)/binomial(k+7,k)',
        objective_definition='f(0)+C; multiply by pi^8*(s/4)^8/Gamma(9)',
        elapsed_seconds=time.monotonic()-start)
    if fit.success:
        expression=float(fit.fun*np.pi**8*(s/4)**8/math.gamma(9))
        result.update(objective=float(fit.fun),density_expression=expression,
            C=float(fit.x[-1]),radial=fit.x[:-1].tolist(),
            max_scaled_residual=float(np.max(A@fit.x)),normalization_residual=float(eq@fit.x-1),
            comparator=COMPARATOR,gap_to_comparator=expression-COMPARATOR)
        np.savez_compressed(D/(name+'_duals.npz'),pair=x,
            scaled_marginals=fit.ineqlin.marginals,scales=scale)
    write(name+'.json',result)
    print('LP',name,{k:v for k,v in result.items() if k not in ['radial','source_nesting']},flush=True)
    return result

def main():
    start=time.monotonic();source=json.loads((I/'radial_result.json').read_text())
    sa,P,H=audit(source,'source63')
    rows=shared_rows(source,P,H)
    control=solve(rows,63,2.,'matched63',source)
    candidate=solve(rows,79,2.,'enriched79',source)
    audits={}
    if candidate['success']:
        audits['enriched79'],_,_=audit(candidate,'enriched79')
    if control['success']:
        audits['matched63'],_,_=audit(control,'matched63')
    # No extra degree-79 search: the two inexpensive scale tests are only a fallback.
    promising=(candidate['success'] and candidate['density_expression']<COMPARATOR-1e-7)
    scales=[]
    if not promising:
        for s in [1.95,2.05]:
            # Preserve pair locations as v=X/(pi*s), clipped to the finite pair annulus.
            x=np.r_[rows['pair']*s/2,np.pi*s,40.]
            x=np.unique(x[(x>=np.pi*s)&(x<=40)])
            sr=dict(rows,pair=x)
            c=solve(sr,63,s,'scale'+str(s).replace('.','_'))
            scales.append(c)
            if c['success'] and c['density_expression']<COMPARATOR:
                audit(c,'scale'+str(s).replace('.','_'))
    write('run_status.json',dict(status='completed',certified_bound=False,
        angular_refits=0,degree79_solves=1,degree63_control_solves=1,scale_solves=len(scales),
        promising_numerical_headroom=promising,elapsed_seconds=time.monotonic()-start))

if __name__=='__main__':main()
