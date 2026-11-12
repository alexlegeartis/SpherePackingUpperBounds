"""Fast exact checks of the useful seed results and rejected proposals."""
import json
import pathlib
import numpy as np
import sympy as sp

root=pathlib.Path(__file__).resolve().parent;p=root/'dim48mod3'
f=np.load(p/'flexes.npz');H,D,V=f['H'],f['Ds'],f['V']
assert np.all(np.einsum('vi,ij,vj->v',V,H,V)==8)
for direction in D:
    assert np.all(np.einsum('vi,ij,vj->v',V,direction,V)==0)
P=np.load(p/'P.npy');coeff=np.array([-2,4,6,4,2,28],np.int64)
assert np.array_equal(618*H+np.einsum('i,ijk->jk',coeff,D),824*P)
assert sp.Matrix(P.tolist()).det()==1
assert np.all(P==P.T) and np.all(P.diagonal()%2==0)
assert sp.Matrix(H.tolist()).det()==3**24
assert np.all(np.einsum('vi,ij,vj->v',V,P,V)==6)
for k,expected in enumerate([sp.Rational(6),sp.Rational(13,2),sp.Rational(5)]):
    d=np.load(root/'CTturyn1'/f'candidate{k}.npz')
    report=json.loads((root/'CTturyn1'/f'candidate{k}_report.json').read_text())
    v=np.array(report['witness'],np.int64)
    assert sp.Rational(int(v@d['Q']@v),int(d['den']))==expected
rejections=json.loads((root/'cluster4_rejections.json').read_text())
for rejected in rejections:
    G=np.array(rejected['G'],np.int64);v=np.array(rejected['rejection_exact_vector'],np.int64)
    assert int(v@G@v)==rejected['rejection_squared_norm']==2
print('PASS:12485 exact norm8 contacts preserve six integer directions; endpoint is even unimodular with exhibited norm6 vectors; all three CT1 proposals and the six-point anticode proposal have exact rejecting witnesses. No minimum or new density claim.')
