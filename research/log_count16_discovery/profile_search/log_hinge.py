"""One nested, dual-priced log-hinge profile experiment. DISCOVERY ONLY.

Retains the seven original v-hinges exactly and adds seven log-hinges.
All cyclic expressions are fully signed. No full proof or density bound is
claimed. Run with: OPENBLAS_NUM_THREADS=1 python log_hinge.py
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog, minimize
from scipy.special import eval_genlaguerre, gammaln

import experiment as base

D = Path(__file__).resolve().parent
V_KNOTS = np.array([1.02, 1.05, 1.1, 1.2, 1.4, 2.0])
LOG_KNOTS = np.array([1.2, 1.4, 1.5, 1.6, 1.65, 1.75])
ALL_KNOTS = np.unique(np.r_[V_KNOTS, LOG_KNOTS])
base.KNOTS = ALL_KNOTS  # Only used by inherited grid, never our basis.
TAIL_SLOPE_FLOOR = 7.0
MARGIN = 2e-7
COMPARATOR = .02492121


def write(name, obj):
    (D/name).write_text(json.dumps(obj, indent=2) + '\n')


class LogHingeExperiment(base.Experiment):
    @staticmethod
    def basis(v):
        v = np.asarray(v)
        return np.concatenate([(v-1)[...,None],
            np.maximum(v[...,None]-V_KNOTS, 0),
            np.log(v)[...,None],
            np.maximum(np.log(v[...,None]/LOG_KNOTS), 0)], axis=-1)

    def profile_info(self, p):
        v = np.unique(np.r_[base.VPRICE, LOG_KNOTS, 1.61532400564816,
                             1.7298028904199105, 1.7333803555690277])
        phi = self.basis(v)@p
        return dict(parameters=p.tolist(), v_knots=V_KNOTS.tolist(),
            log_knots=LOG_KNOTS.tolist(), final_linear_slope=float(p[:7].sum()),
            total_log_slope=float(p[7:].sum()),
            profile_rows=[dict(v=float(x), phi=float(y), h=float(np.exp(-y)),
                              allowance=float((2*np.exp(-y)-np.exp(-2*y))/self.B))
                          for x,y in zip(v,phi)])


def pool_rows(e):
    rows = e.grid(dense=True)
    tt, uu = [], []
    for filename in ['dim16_cyclic_audit.json','dim16_dual_refit_cyclic_audit.json']:
        for c in json.loads((D/'frozen_control'/filename).read_text())['outcomes']:
            tt.append(c['worst_tu'][0]); uu.append(c['worst_tu'][1])
    # Explicit u^2 fibers at every old/new knot, with their exact moving endpoint.
    for k in ALL_KNOTS:
        u = np.sqrt(k)
        t = np.linspace(-1,1/(2*u),1601)
        tt.extend(t); uu.extend(np.full(len(t),u))
    # The sharp isosceles obstacle and each knot side receive local mesh coverage.
    local_t = np.unique(np.r_[np.linspace(.15,.22,7001),
        *(1-ALL_KNOTS/2 + h for h in [-1e-7,0,1e-7])])
    local_t = local_t[(local_t>=-1)&(local_t<=.5)]
    for u in [1.,1.+1e-7,1.+1e-5,1.+1e-3]:
        tt.extend(local_t); uu.extend(np.full(len(local_t),u))
    extra = e.rows(np.array(tt),np.array(uu))
    return {k:np.concatenate([rows[k],extra[k]],axis=0) for k in rows}


def sampled_audit(e,p):
    rows=e.grid(dense=True)
    values=e.expression(rows,p)
    # Independent randomized separation rows, uniform in u and log u.
    rng=np.random.default_rng(16061111)
    u=np.r_[1+3*rng.random(120000), np.exp(np.log(32)*rng.random(120000))]
    t=-1+(1+1/(2*u))*rng.random(len(u))
    rr=e.rows(t,u); vv=e.expression(rr,p)
    candidates=[]
    for r,v in [(rows,values),(rr,vv)]:
        for i in np.argsort(v)[-35:]:
            candidates.append((float(v[i]),float(r['t'][i]),float(r['u'][i])))
    # Separate on every knot fiber as well as in the two-dimensional domain.
    worst=max(candidates)
    local=[]
    for _,t,u in sorted(candidates,reverse=True)[:50]:
        def fun(z):
            u,v=z
            t=-1+v*(1+1/(2*u))
            return -float(e.expression(e.rows(np.array([t]),np.array([u])),p)[0])
        fit=minimize(fun,[u,(t+1)/(1+1/(2*u))],method='Nelder-Mead',
            bounds=[(1,32),(0,1)],options={'maxiter':250,'xatol':2e-12,'fatol':2e-12})
        z=fit.x; t=-1+z[1]*(1+1/(2*z[0]))
        local.append((-float(fit.fun),float(t),float(z[0])))
    for knot in ALL_KNOTS:
        for kind in ['u2','w2']:
            def fv(z):
                a=float(np.asarray(z).reshape(-1)[0])
                if kind=='u2':
                    u=np.sqrt(knot); t=-1+a*(1+1/(2*u))
                else:
                    u=1+a*(np.sqrt(knot)-1); t=(1+u*u-knot)/(2*u)
                return -float(e.expression(e.rows(np.array([t]),np.array([u])),p)[0])
            for start in [0.,.25,.5,.75,1.]:
                fit=minimize(fv,[start],method='Nelder-Mead',bounds=[(0,1)],
                    options={'maxiter':120,'xatol':2e-12,'fatol':2e-12})
                a=float(fit.x[0])
                if kind=='u2':
                    u=np.sqrt(knot); t=-1+a*(1+1/(2*u))
                else:
                    u=1+a*(np.sqrt(knot)-1); t=(1+u*u-knot)/(2*u)
                local.append((-float(fit.fun),float(t),float(u)))
    worst=max([worst]+local)
    return dict(certified=False,dense_rows=len(values),random_rows=len(vv),
        knot_fibers=['u2','w2'],all_knots=ALL_KNOTS.tolist(),
        dense_maximum=float(values.max()),random_maximum=float(vv.max()),
        local_runs=len(local),separated_maximum=worst[0],worst_tu=list(worst[1:]),
        sampled_feasible=worst[0]<=1e-8)


def radial_once(e,p):
    rr=np.load(D/'frozen_control/dim16_radial_rows.npz')
    x,y,z=rr['pair'],rr['fourier'],rr['sign']
    n,d,s=16,63,2.
    norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1))
    sign=(-1.)**np.arange(d+1)
    def lag(X):
        return np.array([eval_genlaguerre(k,n/2-1,2*X)/norm[k] for k in range(d+1)]).T
    h=np.exp(-e.basis(x/(np.pi*s))@p)
    g=np.exp(x)*(2*h-h*h)/e.B
    A=np.vstack([np.c_[lag(x),-g],np.c_[lag(z),np.zeros(len(z))],
                 np.c_[-lag(y)*sign,np.zeros(len(y))]])
    scale=np.maximum(np.abs(A).max(axis=1),1e-200); A/=scale[:,None]
    eq=np.r_[sign,0.]
    fit=linprog(np.ones(d+2),A_ub=A,b_ub=np.full(len(A),-1e-9),
        A_eq=eq[None,:],b_eq=[1.],bounds=[(None,None)]*(d+1)+[(0,None)],method='highs',
        options={'time_limit':120,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    result=dict(certified=False,index='log_hinge',success=bool(fit.success),message=fit.message,
        degree=d,s=s,cut_X=40.,row_count=len(A),candidate_lp_count=1,
        same_rows_as_frozen_control=True,
        row_file_sha256=hashlib.sha256((D/'frozen_control/dim16_radial_rows.npz').read_bytes()).hexdigest())
    if fit.success:
        expression=float(fit.fun*np.pi**8*(s/4)**8/math.gamma(9))
        result.update(objective=float(fit.fun),density_expression=expression,C=float(fit.x[-1]),
            radial=fit.x[:-1].tolist(),max_scaled_residual=float(np.max(A@fit.x)),
            normalization_residual=float(eq@fit.x-1),comparator=COMPARATOR,
            gap_to_comparator=expression-COMPARATOR)
        np.savez_compressed(D/'candidate_pair_dual_pricing.npz',pair=x,
            scaled_marginals=fit.ineqlin.marginals[:len(x)],fixed_row_scale=scale[:len(x)])
    return result


def main():
    started=time.monotonic()
    e=LogHingeExperiment(16)
    source=json.loads((D/'frozen_control/dim16_optimization.json').read_text())['candidates'][0]
    control=json.loads((D/'frozen_control/dim16_dual_refit.json').read_text())['candidate']
    p0=np.r_[control['parameters'],np.zeros(7)]
    prices=np.load(D/'frozen_control/dim16_pair_dual_pricing.npz')
    gradient=prices['tangent_weight']@e.basis(prices['pair']/(2*np.pi))
    raw_gradient=gradient.copy(); gradient/=np.linalg.norm(gradient)
    rows=pool_rows(e)
    # Validate exact numerical nesting separately from the changed grid.
    v=np.r_[1.,ALL_KNOTS,np.linspace(1,8,1001)]
    oldbasis=np.c_[v-1,np.maximum(v[:,None]-V_KNOTS,0)]
    nesting_error=float(np.max(np.abs(e.basis(v)@p0-oldbasis@np.array(control['parameters']))))
    scale=np.maximum(1e-3,np.abs(rows['q0'])+np.abs(rows['q1'])+np.abs(rows['q2']))
    cv=e.expression(rows,p0)
    print('CONTROL',dict(rows=len(cv),maximum=float(cv.max()),nesting_error=nesting_error),flush=True)
    con={'type':'ineq','fun':lambda p:(-e.expression(rows,p)-MARGIN)/scale,
         'jac':lambda p:-e.expression(rows,p,True)[1]/scale[:,None]}
    tail={'type':'ineq','fun':lambda p: p[:7].sum()-TAIL_SLOPE_FLOOR,
          'jac':lambda p:np.r_[np.ones(7),np.zeros(7)]}
    fit=minimize(lambda p:float(gradient@p),p0,jac=lambda p:gradient,method='SLSQP',
        bounds=[(0,100)]*14,constraints=[con,tail],
        options={'maxiter':120,'ftol':2e-12,'disp':True})
    values=e.expression(rows,fit.x); iw=int(np.argmax(values))
    out=dict(certified=False,dimension=16,source_sha256=e.sha256,fit_count=1,
        pricing_source_index=0,start_and_control_index='dual_refit',
        raw_tangent_gradient=raw_gradient.tolist(),normalized_tangent_gradient=gradient.tolist(),
        training_rows=len(values),full_signed_cyclic_rows=True,margin=MARGIN,
        tail_slope_floor=TAIL_SLOPE_FLOOR,nesting_error=nesting_error,
        success=bool(fit.success),message=fit.message,iterations=int(fit.nit),
        objective=float(fit.fun),control_objective=float(gradient@p0),
        train_maximum=float(values[iw]),train_worst_tu=[float(rows['t'][iw]),float(rows['u'][iw])],
        elapsed_seconds=time.monotonic()-started,**e.profile_info(fit.x))
    write('fit.json',out)
    print('FIT',{k:v for k,v in out.items() if k not in ['profile_rows','raw_tangent_gradient','normalized_tangent_gradient']},flush=True)
    if fit.success and values.max()<=1e-8:
        audit=sampled_audit(e,fit.x); write('cyclic_audit.json',audit)
        print('AUDIT',audit,flush=True)
        if audit['sampled_feasible']:
            radial=radial_once(e,fit.x); write('radial_result.json',radial)
            print('RADIAL',{k:v for k,v in radial.items() if k!='radial'},flush=True)
    write('run_status.json',dict(status='completed',certified=False,elapsed_seconds=time.monotonic()-started,
        full_proof_attempted=False,fit_count=1,radial_candidate_runs=int((D/'radial_result.json').exists())))


if __name__=='__main__':
    main()
