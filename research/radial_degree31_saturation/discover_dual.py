"""Zero-margin LP selects rows only; this script is not a certificate."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from flint import arb, ctx, fmpq as Q
from radial_model import Model,sample_specs,D

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
ap=argparse.ArgumentParser()
ap.add_argument('--dimension',type=int,required=True)
ap.add_argument('--B')
ap.add_argument('--tag',default='fixed')
args=ap.parse_args()
ctx.prec=640
n=args.dimension
B=args.B or json.loads((D/f'inputs/angular{n}.json').read_text())['B_exact']
model=Model(n,B)
specs=sample_specs(n)
A=[];a0=[]
for spec in specs:
    _,_,row=model.row(spec,False)
    largest=max(float(abs(x).upper()) for x in row)
    exponent=math.frexp(largest)[1] if largest else 0
    multiplier=Q(1,2**exponent) if exponent>=0 else Q(2**(-exponent))
    spec['multiplier_exact']=str(multiplier)
    transformed,constant,_=model.row(spec)
    A.append([float(x.mid()) for x in transformed]);a0.append(float(constant.mid()))
A=np.array(A);a0=np.array(a0)
c=np.array(list(map(float,model.c)))
res=linprog(c,A_ub=A,b_ub=-a0,bounds=[(0,None)]+[(None,None)]*31+[(0,None)],method='highs',
            options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
assert res.success,res.message
y=-res.ineqlin.marginals
active=np.flatnonzero(y>0)
out={'dimension':n,'B_exact':str(Q(B)),'tag':args.tag,'objective_float':1+float(res.fun),
     'density_float':(1+float(res.fun))*math.pi**(n/2)*1.5**(n/2)/math.gamma(n/2+1),
     'active_count':len(active),'minimum_active_y':float(min(y[active])),
     'stationarity_residual_float':float(np.max(abs(A.T@y+c))),
     'nonnegative_variable_reduced_costs_float':res.lower.marginals[[0,-1]].tolist(),
     'primal_t_C_float':res.x[[0,-1]].tolist(),'row_count':len(specs),
     'rows':[specs[i] for i in active],'y_float':y[active].tolist(),
     'primal_float':res.x.tolist(),'certified':False}
dest=D/f'dual{n}_{args.tag}_selected.json'
dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['rows','y_float','primal_float']},indent=2),flush=True)
