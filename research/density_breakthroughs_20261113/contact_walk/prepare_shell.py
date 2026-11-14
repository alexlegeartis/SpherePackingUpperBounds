"""Prepare a floating discovery search from exact Gram data; never a proof."""
from pathlib import Path
import numpy as np
import json
P=Path(__file__).resolve().parent
S=P.parents[2]/'constructions/dim36'
Q=np.loadtxt(S/'Q_reduced.txt',dtype=np.int64)
U=np.loadtxt(S/'data/U.txt',dtype=np.int64)
np.savetxt(P/'Q.txt',Q,fmt='%d')
np.savetxt(P/'R.txt',np.linalg.cholesky(Q).T,fmt='%.17g')
np.save(P/'U.npy',U)
print(json.dumps({'dimension':36,'threshold':160,'scope':'floating discovery with exact leaf norms'}))
