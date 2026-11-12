#!/usr/bin/env python3
"""Reconstruct explicit GF(4) outer codes for mixed hexagonal packings.

GF(4) is encoded as 0,1,w,w+1 = 0,1,2,3, with w^2=w+1.
Sources are archived HTML and the existing, checked 36D provenance.
"""
import json
import re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'constructions/dim36/provenance/ct36_symmetry_break_20261112'
MUL = np.array([[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]], dtype=np.uint8)

def shorten(g, positions):
    g = g.copy()
    r = 0
    for c in positions:
        pivot = next(i for i in range(r,len(g)) if g[i,c])
        g[[pivot,r]] = g[[r,pivot]]
        inv = next(i for i in range(1,4) if MUL[i,g[r,c]] == 1)
        g[r] = MUL[inv,g[r]]
        for i in range(len(g)):
            if i != r and g[i,c]:
                g[i] ^= MUL[g[i,c],g[r]]
        r += 1
    assert not g[r:,positions].any()
    return np.delete(g[r:], positions, axis=1)

def write(name, g):
    np.save(HERE / f'{name}_G.npy',g)
    packed = []
    for row in g:
        for a in (1,2):
            packed.append(sum(int(x) << (2*i) for i,x in enumerate(MUL[a,row])))
    (HERE/f'{name}.input').write_text(f'{g.shape[1]} {2*len(g)}\n'+' '.join(map(str,packed))+'\n')

def main():
    g20 = np.array(json.loads((OLD/'outer18.json').read_text())['G20'],dtype=np.uint8)
    write('n20k13',g20)
    write('n19k12',shorten(g20,[19]))
    write('n18k11',shorten(g20,[18,19]))
    html = (HERE/'sources/q4_n17_k9.html').read_text()
    rows = re.findall(r'^\[ ([01w^2, ]+) \]',html,re.M)
    assert len(rows)==9
    enc = {'0':0,'1':1,'w':2,'w^2':3}
    g18 = np.array([[enc[x.strip()] for x in row.split(',')] for row in rows],dtype=np.uint8)
    assert g18.shape==(9,18)
    write('n17k9',g18[:,:17])

if __name__=='__main__':
    main()
