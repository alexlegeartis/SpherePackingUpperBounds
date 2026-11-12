import numpy as np, json
from pathlib import Path
P=Path(__file__).resolve().parent
A=np.load(P.parent/'ct36_metric_screen_20261112/parts.npy')
N=A[0]+2*A[1]+5*A[2]
g=int(np.gcd.reduce(N.flatten())); print("gcd",g)
H=N/36864;np.save(P/"N.npy",N//g)
np.save(P/'H.npy',H)
np.savetxt(P/'R.txt',np.linalg.cholesky(H).T)
print('integral',np.diag(H),np.linalg.slogdet(H),flush=True)
