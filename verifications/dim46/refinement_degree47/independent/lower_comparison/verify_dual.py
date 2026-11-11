"""Validated dual LOWER certificate for a fixed finite-dimensional radial family.

Selected rows are reevaluated with true pi and exponential functions. Arb's
verified solve encloses the exact solution of the true active-row system;
floating residuals and dyadic approximations are never used as certificates.
"""
import argparse
import json
import math
from pathlib import Path
from flint import arb, arb_mat, fmpq as Q, ctx
from radial_model import Model,D

if not __debug__:
    raise RuntimeError('Proof checks require assertions; do not use Python -O')
ap=argparse.ArgumentParser()
ap.add_argument('selected')
ap.add_argument('--bits',type=int,default=1024)
args=ap.parse_args()
ctx.prec=args.bits
path=Path(args.selected)
record=json.loads(path.read_text())
model=Model(record['dimension'],record['B_exact'])
assert len(record['rows'])==33, 'A 33-row square system is required for this first verifier'
rows=[];constants=[]
for spec in record['rows']:
    row,A0,_=model.row(spec)
    rows.append(row);constants.append(A0)
A=arb_mat([[rows[i][j] for i in range(33)] for j in range(33)])
rhs=arb_mat([[-arb(c)] for c in model.c])
determinant=A.det()
assert determinant.is_finite() and not determinant.contains(0)
y=A.solve(rhs,algorithm='precond')
weights=[y[i,0] for i in range(33)]
assert all(w.is_finite() for w in weights)
assert all(w>0 for w in weights), [str(w) for w in weights if not w>0]
lower=1+sum((w*a for w,a in zip(weights,constants)),arb(0))
assert lower>0
n=record['dimension']
volume=(arb.pi()*arb(Q(3,2)))**arb(Q(n,2))/arb(Q(n+2,2)).gamma()
density=lower*volume
digits=16
safe=Q(math.floor(float(density.lower())*10**digits)-1,10**digits)
assert safe>0 and density>arb(safe)
result={'pass':True,'dimension':n,'B_exact':record['B_exact'],'precision_bits':ctx.prec,
        'family':'t>=0 times normalized degree111 source plus free normalized Laguerre degrees0..31; C=BM>=0; H(0)=1; Fourier H>=0; physical P<=0 for X>=7pi; pair coupling on squared radii[6,7]',
        'necessary_rows_only':True,'genuine_pi_and_exponential_rows':True,
        'selected_row_count':len(weights),'all_true_dual_weights_strictly_positive':True,
        'all_32_Laguerre_Fourier_eigenidentities_exact':True,
        'all_pair_rows_at_contact':all(s['kind']!='pair' or Q(s['u'])==6 for s in record['rows']),
        'active_matrix_determinant_interval':str(determinant),
        'minimum_weight_interval':str(min(weights,key=lambda w:float(w.mid()))),
        'objective_lower_interval':str(lower),'density_lower_interval':str(density),
        'safe_density_lower_exact':str(safe),'safe_density_lower_decimal':f'{float(safe):.16f}',
        'dual_definition':'unique exact solution Atilde_active^T y=-c; verified nonsingular Arb solve encloses true weights',
        'nonnegative_variable_reduced_costs_exact':'0 for t and C',
        'weights_intervals':list(map(str,weights)),
        'scope':'Lower bound for this fixed finite-dimensional radial family at the specified B. Not an unrestricted packing-density lower or upper bound.'}
dest=path.with_name(path.stem.replace('_selected','_verified')+'.json')
dest.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='weights_intervals'},indent=2),flush=True)
