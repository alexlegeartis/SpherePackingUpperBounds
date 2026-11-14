"""Independent exact replay of the fixed [17,9,7] rank-two obstruction.

No use is made of the C++ enumerator, its supports, or the main verify.py when
constructing the code.  The archived support list is compared afterward.
"""
from pathlib import Path
from itertools import product
from fractions import Fraction
import hashlib, json
import numpy as np
P=Path(__file__).resolve().parent
M=P.parent/'mixed_codes';g=np.load(M/'n17k9_G.npy')
assert g.shape==(9,17) and np.all((g>=0)&(g<4))
def mul(a,b):
    z=0
    for k in range(2):
        if (a>>k)&1:z^=int(b)<<k
    if z&4:z^=7  # Reduction by x^2+x+1.
    return z
C=np.zeros((1,17),dtype=np.uint8)
for row in g:
    C=np.concatenate([C ^ np.array([mul(a,b)for b in row],dtype=np.uint8) for a in range(4)])
w=np.count_nonzero(C,axis=1)
assert len(C)==4**9 and (w==0).sum()==1 and min(w[w>0])==7
word_masks=np.sum((C!=0)*(1<<np.arange(17,dtype=np.int64)),axis=1)
supports=np.unique(word_masks[w==7]);assert len(supports)==408
stored=np.loadtxt(M/'n17k9.supports',dtype=np.int64,skiprows=1)
assert np.array_equal(supports,np.sort(stored))
B=((supports[:,None]>>np.arange(17))&1).astype(np.int64)
incidence=B.T@B
assert np.all(np.diag(incidence)==168)
assert np.all(incidence[np.triu_indices(17,1)]==63)
colors=[(1,0),(0,1),(1,1),(1,2)]
admissible=[]
for counts in product(range(8),repeat=4):
    if sum(counts)!=7:continue
    sums={(0,0)}
    for color,count in zip(colors,counts):
        for _ in range(count):
            sums={((x+sign*color[0])%3,(y+sign*color[1])%3) for x,y in sums for sign in [-1,1]}
    if (0,0) not in sums:admissible.append(counts)
assert len(admissible)==12 and {tuple(sorted(c))for c in admissible}=={(0,0,1,6)}
required=Fraction(408*15,63);assert required.denominator!=1
result={
 'status':'independently verified exact obstruction',
 'scope':'this fixed GF4 outer code; all parity matrices over F3 of rank at most two with no zero column',
 'generator_sha256':hashlib.sha256((M/'n17k9_G.npy').read_bytes()).hexdigest(),
 'codewords':len(C),'weight_enumerator':np.bincount(w,minlength=18).tolist(),
 'weight7_words':int((w==7).sum()),'weight7_supports':len(supports),
 'design_v_k_lambda':[17,7,63],'replication_number':168,
 'admissible_projective_color_counts':[list(c)for c in admissible],
 'monochromatic_pairs_per_block':15,'required_total_monochromatic_coordinate_pairs':str(required),
 'argument':'Every block would have color counts (6,1), hence 15 monochromatic pairs. Each coordinate pair occurs in 63 blocks. Therefore the integer number of monochromatic coordinate pairs would be 408*15/63=680/7, a contradiction.',
 'zero_column_qualification':'The combinatorial statement explicitly excludes zero columns. In the intended mixed construction, a zero parity column admits a ternary weight-one word of local squared cost4, preventing the required minimum8.'}
(P/'audit_n17.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
