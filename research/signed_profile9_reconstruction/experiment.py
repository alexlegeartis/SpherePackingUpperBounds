"""Fixed-angular, kinked convex log radial profile: numerical discovery only.

Every cyclic row retains all three signed Q contributions. No numerical
diagnostic in this script is a global certificate or a packing upper bound.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog, minimize
from scipy.special import eval_genlaguerre, gammaln

D = Path(__file__).resolve().parent
KNOTS = np.array([1.02, 1.05, 1.1, 1.2, 1.4, 2.0])
VPRICE = np.array([1.005, 1.01, 1.02, 1.035, 1.05, 1.075, 1.1,
                   1.15, 1.2, 1.3, 1.4, 1.6, 2.0])


def write(name, obj):
    (D / name).write_text(json.dumps(obj, indent=2) + '\n')


class Experiment:
    def __init__(self, n=16):
        self.n = n
        path = D / 'inputs' / f'angular_dim{n}.json'
        data = path.read_bytes()
        self.sha256 = hashlib.sha256(data).hexdigest()
        j = json.loads(data)
        self.c = np.array([float(Fraction(z)) for z in j['gegenbauer_coefficients']])
        self.B = float(Fraction(j['B_exact']))
        self.a = float(Fraction(j['a_exact']))
        self.b = float(Fraction(j['b_exact']))

    def q(self, t):
        t = np.asarray(t)
        p, r = np.ones_like(t), t.copy()
        out = self.c[0] * p + self.c[1] * r
        for k in range(1, len(self.c) - 1):
            p, r = r, ((2*k+self.n-2)*t*r-k*p)/(k+self.n-2)
            out += self.c[k+1] * r
        return out

    @staticmethod
    def basis(v):
        v = np.asarray(v)
        return np.concatenate([(v-1)[..., None],
                               np.maximum(v[..., None]-KNOTS, 0)], axis=-1)

    def control_phi(self, v):
        return self.a*(v-1) + self.b*(v-1)**2

    def rows(self, t, u):
        t, u = np.asarray(t), np.asarray(u)
        w2 = 1+u*u-2*u*t
        w = np.sqrt(w2)
        return dict(t=t, u=u, y=u*u, z=w2, q0=self.q(t),
                    q1=self.q((1-u*t)/w), q2=self.q((u-t)/w),
                    diff=self.basis(w2)-self.basis(u*u), end=self.basis(w2))

    def grid(self, dense=False):
        nt = 801 if dense else 321
        us = np.unique(np.r_[np.linspace(1, 1.08, 121 if dense else 41),
                              np.linspace(1.08, 1.5, 101 if dense else 41),
                              np.linspace(1.5, 3, 101 if dense else 31),
                              np.geomspace(3, 32, 91 if dense else 25),
                              np.sqrt(KNOTS)])
        t0 = np.linspace(-1, .5, nt)
        t, u = np.meshgrid(t0, us)
        m = 2*t*u <= 1+1e-14
        t, u = t[m], u[m]
        # Include sorted-equal-largest-side boundary and kink surfaces.
        tx = np.linspace(-1, .5, 1501 if dense else 501)
        ux = np.r_[np.ones_like(tx), us, us]
        tx = np.r_[tx, 1/(2*us), -np.ones_like(us)]
        for z in KNOTS:
            uk = np.linspace(1, np.sqrt(z), 301 if dense else 101)
            tk = (1+uk*uk-z)/(2*uk)
            tx = np.r_[tx, tk]
            ux = np.r_[ux, uk]
        return self.rows(np.r_[t, tx], np.r_[u, ux])

    def expression(self, rows, p=None, jac=False):
        if p is None:
            e1 = np.exp(-(self.control_phi(rows['z'])-self.control_phi(rows['y'])))
            e2 = np.exp(-self.control_phi(rows['z']))
        else:
            e1 = np.exp(-rows['diff']@p)
            e2 = np.exp(-rows['end']@p)
        c1, c2 = e1*rows['q1'], e2*rows['q2']
        value = rows['q0']+c1+c2
        if jac:
            return value, -c1[:, None]*rows['diff']-c2[:, None]*rows['end']
        return value

    def profile_info(self, p):
        phi = self.basis(VPRICE)@p
        h = np.exp(-phi)
        hc = np.exp(-self.control_phi(VPRICE))
        allow, ac = 2*h-h*h, 2*hc-hc*hc
        return dict(parameters=p.tolist(), knots=KNOTS.tolist(),
                    cumulative_slopes=np.cumsum(p).tolist(),
                    allowance_rows=[dict(v=float(v), phi=float(f), control_phi=float(fc),
                        allowance=float(a/self.B), control_allowance=float(b/self.B),
                        relative_gain=float(a/b-1)) for v,f,fc,a,b in zip(VPRICE,phi,self.control_phi(VPRICE),allow,ac)])

    def optimize(self, trials=3):
        rows = self.grid()
        base = self.expression(rows)
        print('GRID', len(base), 'control maximum', base.max(), flush=True)
        edges = np.r_[1., KNOTS]
        slopes = self.a+self.b*((edges[:-1]-1)+(edges[1:]-1))
        # The final interval uses a tangent at v=2; all increments nonnegative.
        slopes = np.r_[slopes, self.a+2*self.b*(KNOTS[-1]-1)]
        start = np.r_[slopes[0], np.diff(slopes)]
        scale = np.maximum(1e-3, np.abs(rows['q0'])+np.abs(rows['q1'])+np.abs(rows['q2']))
        bp = self.basis(VPRICE)
        hctrl = np.exp(-self.control_phi(VPRICE))
        actrl = 2*hctrl-hctrl*hctrl
        results = []
        weights = [np.exp(-((VPRICE-1)/.16)**2),
                   np.exp(-((VPRICE-1)/.08)**2),
                   np.exp(-((VPRICE-1)/.30)**2)]
        for k, ww in enumerate(weights[:trials]):
            ww /= ww.sum()
            def fun(p):
                h = np.exp(-bp@p)
                return -float(np.sum(ww*(2*h-h*h)/actrl))
            def jac(p):
                h = np.exp(-bp@p)
                return np.sum((ww*2*h*(1-h)/actrl)[:, None]*bp, axis=0)
            con = {'type':'ineq',
                   'fun':lambda p: (-self.expression(rows,p)-1e-7)/scale,
                   'jac':lambda p: -self.expression(rows,p,True)[1]/scale[:,None]}
            fit = minimize(fun, start, jac=jac, method='SLSQP', bounds=[(0, 100)]*7,
                           constraints=[con], options={'maxiter':120,'ftol':2e-12,'disp':False})
            value = self.expression(rows,fit.x)
            j = int(np.argmax(value))
            out = dict(index=k, success=bool(fit.success), message=fit.message,
                       iterations=int(fit.nit), objective=float(fit.fun),
                       train_maximum=float(value[j]), train_worst=dict(t=float(rows['t'][j]),u=float(rows['u'][j])),
                       **self.profile_info(fit.x))
            results.append(out)
            write(f'dim{self.n}_optimization.json', dict(dimension=self.n, certified=False,
                    source_sha256=self.sha256,B=self.B,training_row_count=len(base),
                    control_train_maximum=float(base.max()), candidates=results))
            print('FIT', k, {z:out[z] for z in ['success','message','iterations','objective','train_maximum','parameters']},flush=True)
        return results

    def audit(self, candidates):
        rows = self.grid(dense=True)
        outcomes=[]
        for cc in [None]+candidates:
            p = None if cc is None else np.array(cc['parameters'])
            vals = self.expression(rows,p)
            ids = np.argsort(vals)[-30:]
            worst = float(vals.max())
            point = [float(rows['t'][ids[-1]]),float(rows['u'][ids[-1]])]
            # One bounded continuous separation pass, including kink sides.
            for j in ids:
                def f(z):
                    u,v=z
                    t=-1+v*(1+1/(2*u))
                    return -float(self.expression(self.rows(np.array([t]),np.array([u])),p)[0])
                u0=rows['u'][j]
                v0=(rows['t'][j]+1)/(1+1/(2*u0))
                opt=minimize(f,[u0,v0],method='Nelder-Mead',bounds=[(1,32),(0,1)],
                             options={'maxiter':150,'xatol':1e-11,'fatol':1e-11})
                if -opt.fun>worst:
                    worst=float(-opt.fun);point=[float(-1+opt.x[1]*(1+1/(2*opt.x[0]))),float(opt.x[0])]
            outcome=dict(index='quadratic' if cc is None else cc['index'],dense_rows=len(vals),
                         dense_maximum=float(vals.max()),separated_maximum=worst,worst_tu=point,
                         sampled_feasible=worst<=1e-8)
            outcomes.append(outcome)
            print('AUDIT',outcome,flush=True)
        write(f'dim{self.n}_cyclic_audit.json',dict(certified=False,outcomes=outcomes))
        return outcomes

    def radial(self, candidates, degree=63, filename=None):
        n,d,s,cut = self.n,degree,self.n/8,40.
        norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1))
        sign=(-1.)**np.arange(d+1)
        def lag(X):
            return np.array([eval_genlaguerre(k,n/2-1,2*X)/norm[k] for k in range(d+1)]).T
        y=np.r_[np.linspace(0,4*d+n+20,6001),np.geomspace(4*d+n+20,2500,801)]
        x=np.unique(np.r_[np.linspace(np.pi*s,cut,1501),np.pi*s+np.geomspace(1e-8,2,501)])
        z=np.r_[np.linspace(cut,4*d+n+20,4001),np.geomspace(4*d+n+20,2500,801)]
        lx, ly, lz=lag(x),lag(y),lag(z)
        eq=np.r_[sign,0.]
        output=[]
        np.savez_compressed(D/f'dim{n}_radial_rows.npz',pair=x,fourier=y,sign=z)
        for cc in [None]+candidates:
            phi=self.control_phi(x/(np.pi*s)) if cc is None else self.basis(x/(np.pi*s))@np.array(cc['parameters'])
            h=np.exp(-phi)
            g=np.exp(x)*(2*h-h*h)/self.B
            A=np.vstack([np.c_[lx,-g],np.c_[lz,np.zeros(len(z))],np.c_[-ly*sign,np.zeros(len(y))]])
            scale=np.maximum(np.abs(A).max(axis=1),1e-200)
            A/=scale[:,None]
            fit=linprog(np.ones(d+2),A_ub=A,b_ub=np.full(len(A),-1e-9),A_eq=eq[None,:],b_eq=[1.],
                        bounds=[(None,None)]*(d+1)+[(0,None)],method='highs',
                        options={'time_limit':120,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
            out=dict(index='quadratic' if cc is None else cc['index'],success=bool(fit.success),message=fit.message)
            if fit.success:
                out.update(objective=float(fit.fun),density_expression=float(fit.fun*np.pi**(n/2)*(s/4)**(n/2)/math.gamma(n/2+1)),
                           C=float(fit.x[-1]),radial=fit.x[:-1].tolist(),max_scaled_residual=float(np.max(A@fit.x)),
                           normalization_residual=float(eq@fit.x-1))
            output.append(out)
            print('RADIAL',{k:v for k,v in out.items() if k!='radial'},flush=True)
            write(filename or f'dim{n}_radial_comparison.json',dict(dimension=n,degree=d,s=s,cut_X=cut,B=self.B,
                  certified=False,row_count=len(A),same_pair_fourier_sign_rows=True,outcomes=output))
        return output


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--dimension',type=int,default=16)
    ap.add_argument('--stage',choices=['optimize','audit','radial','all'],default='all')
    ap.add_argument('--trials',type=int,default=3)
    args=ap.parse_args()
    e=Experiment(args.dimension)
    if args.stage in ['optimize','all']:
        candidates=e.optimize(args.trials)
    else:
        candidates=json.loads((D/f'dim{e.n}_optimization.json').read_text())['candidates']
    if args.stage in ['audit','all']:
        audits=e.audit(candidates)
    else:
        audits=json.loads((D/f'dim{e.n}_cyclic_audit.json').read_text())['outcomes'] if args.stage=='radial' else []
    if args.stage in ['radial','all']:
        eligible=[c for c in candidates if any(a['index']==c['index'] and a['sampled_feasible'] for a in audits)]
        e.radial(eligible)
