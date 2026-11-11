"""One conditional U=7803 count test and at most one bounded root-cut refit.

No spherical count proof is assumed completed. No density bound is asserted.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,time,math,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
import radial_core as e

D=Path(__file__).resolve().parent
TAU=F(10000,4997);U=7803;W=7333;S=2.;DEG=63
VOL=np.pi**8*(S/4)**8/math.gamma(9)

def write(name,obj):e.write(name,obj)

def initial_rows():
    z=np.load(D/'inputs/shared_rows.npz');rows={k:z[k].copy() for k in z.files}
    for name in ['source63','matched63','enriched79']:
        a=json.loads((D/'inputs'/f'{name}_audit.json').read_text())
        for gate,key in [('fourier','fourier'),('physical','sign')]:
            rows[key]=np.r_[rows[key],[w['x'] for w in a['gates'][gate]['counterexamples']]]
    rows['pair']=np.r_[rows['pair'],np.pi*np.linspace(2,2.015,4001)]
    return {k:np.unique(v) for k,v in rows.items()}

def solve(rows,name,margin):
    start=time.monotonic()
    x=np.unique(rows['pair']);y=rows['fourier'];z=rows['sign']
    # Explicit rational threshold-side metadata takes precedence over float comparison.
    delta=F(1,2**30)
    special=[('left',TAU-delta),('at',TAU),('right',TAU+delta)]
    active=x<np.pi*float(TAU)
    side_rows=[]
    for side,point in special:
        side_rows.append(dict(side=side,row=len(x),squared_radius_exact=str(point),active=bool(point<TAU)))
        x=np.r_[x,np.pi*float(point)];active=np.r_[active,point<TAU]
    h=np.exp(-e.basis(x/(np.pi*S))@e.PARS)
    g=np.exp(x)*(2*h-h*h)/e.B
    countcol=np.exp(x)*active/W
    sign=(-1.)**np.arange(DEG+1)
    A=np.vstack([np.c_[e.lag(x,DEG),-g,-countcol],
        np.c_[e.lag(z,DEG),np.zeros((len(z),2))],
        np.c_[-e.lag(y,DEG)*sign,np.zeros((len(y),2))]])
    scale=np.maximum(abs(A).max(axis=1),1e-200);A/=scale[:,None]
    eq=np.r_[sign,0.,0.];objective=np.r_[np.ones(DEG+2),U/W]
    b=np.full(len(A),-margin)
    np.savez_compressed(D/(name+'_rows.npz'),pair=x,fourier=y,sign=z,active=active,
        A=A,b=b,scale=scale,equality=eq,objective=objective)
    fit=linprog(objective,A_ub=A,b_ub=b,A_eq=eq[None,:],b_eq=[1.],
        bounds=[(None,None)]*(DEG+1)+[(0,None),(0,None)],method='highs',
        options={'time_limit':120.,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    result=dict(certified_bound=False,count_conditional=True,count_proof_received=False,
        success=bool(fit.success),message=fit.message,dimension=16,degree=DEG,s=S,cut_X=40,
        U=U,threshold_squared_exact=str(TAU),radius_ratio_squared_exact=str(TAU/F(2)),
        half_open_shell='2 <= r^2 < 10000/4997',
        pair_constraint='exp(-X)*P(X) <= C/B*(2*h-h^2) + lambda*1[2 <= X/pi < 10000/4997]',
        objective_definition='sum(radial)+C+7803*lambda; multiplied by pi^8*(2/4)^8/Gamma(9)',
        count_variable='W=7333*lambda; ordinary average-neighbor count column, not a PSD kernel',
        row_counts=dict(pair=len(x),fourier=len(y),physical=len(z)),margin=margin,
        threshold_rows=side_rows,elapsed_seconds=time.monotonic()-start)
    if fit.success:
        expression=float(fit.fun*VOL)
        result.update(objective=float(fit.fun),density_expression=expression,
            comparator=e.COMPARATOR,gap_to_comparator=expression-e.COMPARATOR,
            radial=fit.x[:DEG+1].tolist(),C=float(fit.x[-2]),lambda_count=float(fit.x[-1]/W),
            count_objective_cost=float(U*fit.x[-1]/W),
            max_scaled_residual=float(np.max(A@fit.x)),normalization_residual=float(eq@fit.x-1))
    write(name+'.json',result)
    print('CONDITIONAL LP',name,{k:v for k,v in result.items() if k not in ['radial','threshold_rows']},flush=True)
    return result

def audit(case,name):
    result,P,H=e.audit(case,name)
    export=json.loads((D/(name+'_exact.json')).read_text())
    export.update(lambda_count_exact=str(e.rat(case['lambda_count'])),U_exact=str(U),
        threshold_squared_exact=str(TAU),count_conditional=True)
    write(name+'_exact.json',export)
    return result,P,H

def refine(rows,diagnostic,P,H):
    detail={}
    for gate,key,p,lo in [('fourier','fourier',H,0),('physical','sign',P,40)]:
        critical=e.roots(p.derivative(),lo)+e.roots(p.derivative()-p,lo)
        points=[r['x'] for r in critical]
        rts=[r['x'] for r in diagnostic['gates'][gate]['roots']]
        intervals=[]
        for w in diagnostic['gates'][gate]['counterexamples']:
            left=max([float(lo)]+[x for x in rts if x<w['x']])
            right=min([x for x in rts if x>w['x']]+[max(rts+[2500])+1])
            width=right-left
            # Bounded neighborhood, dense enough to challenge movement of a near-double root.
            l=max(float(lo),left-width);r=right+width
            points.extend(np.linspace(l,r,1501))
            intervals.append(dict(left=l,right=r,rows=1501))
        rows[key]=np.unique(np.r_[rows[key],points])
        detail[gate]=dict(critical_roots=critical,neighborhoods=intervals)
    write('root_cut_refinement.json',dict(detail=detail,new_counts={k:len(v) for k,v in rows.items()},
        scope='One bounded refit only; new critical points plus 1501-point expanded neighborhoods around each bad interval.'))
    return rows

def main():
    started=time.monotonic();rows=initial_rows()
    initial=solve(rows,'conditional_initial',1e-9)
    did_refit=False
    if initial['success']:
        aa,P,H=audit(initial,'conditional_initial')
        needs_refit=any(not v['pass_conservative_root_gate'] for v in aa['gates'].values())
        if needs_refit:
            did_refit=True;rr=refine(rows,aa,P,H)
            final=solve(rr,'conditional_refit',1e-8)
            if final['success']:audit(final,'conditional_refit')
    write('run_status.json',dict(status='completed',certified_bound=False,count_conditional=True,
        initial_solves=1,root_cut_refits=int(did_refit),angular_refits=0,
        elapsed_seconds=time.monotonic()-started))

if __name__=='__main__':main()
