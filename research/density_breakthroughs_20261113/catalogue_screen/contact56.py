#!/usr/bin/env python3
"""Exact eutaxy identity and numerical perfection screen for the 56D seed."""
import json,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh
P=Path(__file__).resolve().parent
d=np.load(P/'L56a_orbit.npz');G=d['G'];V=d['V'];F=d['frame']
assert np.array_equal(F,26970*G)
rng=np.random.default_rng(560013)
W=V[rng.choice(len(V),6000,False)]@np.linalg.cholesky(G)
i,j=np.triu_indices(56);A=W[:,i]*W[:,j]*np.where(i==j,1,2)
t=time.time();e=eigvalsh(A.T@A)
report={'sample_contacts':len(W),'symmetric_variables':len(i),'numerical_rank':int(sum(e>e[-1]*1e-9)),'smallest_eigenvalues':e[:10].tolist(),'largest_eigenvalue':float(e[-1]),'exact_eutaxy_identity':'G (sum_v vv^T) G = 26970 G','identity_checked':True,'scope':'sampled numerical perfection rank; exact positive eutaxy identity','seconds':time.time()-t}
(P/'L56a_contact_screen.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
