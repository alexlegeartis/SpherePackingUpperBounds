"""Explicit primitive 53D section of Leterrier's published 54D Gram.

The parent minimum10 is inherited from Leterrier2024; all source identity,
normal, kernel, determinant and minimum-attainment checks below are exact.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd
import argparse,json,hashlib,subprocess,re
import numpy as np
from flint import fmpz_mat,fmpq_mat
from common import P,finish
if not __debug__:raise RuntimeError('Run without -O.')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()


def audit(source_pdf=None):
    source=P/'sources/leterrier54_gram.txt';sources=json.loads((P/'sources.json').read_text())
    assert sha(source)==sources['Leterrier2024']['matrix_sha256']
    rows=[list(map(int,line.split()))for line in source.read_text().splitlines()]
    assert len(rows)==54 and all(len(v)==54 for v in rows)
    if source_pdf:
        assert sha(source_pdf)==sources['Leterrier2024']['pdf_sha256']
        text=subprocess.check_output(['pdftotext','-layout',str(source_pdf),'-'],text=True)
        text=text.split('Computational proposition 3.3.10.',1)[1].split('Proof. –– 1.',1)[0]
        parsed=[list(map(int,re.findall(r'-?\d+',line.replace('−','-'))))for line in text.splitlines()]
        assert [v for v in parsed if len(v)>20]==rows
    with np.load(P/'leterrier53.npz',allow_pickle=False)as d:
        assert all(np.issubdtype(d[k].dtype,np.integer)for k in d.files)
        assert d['G'].tolist()==rows
        a=[int(x)for x in d['normal_label']];w=[int(x)for x in d['minimum_witness']]
    G=fmpz_mat(rows);assert G==G.transpose()and G.det()==3**25
    assert all(fmpz_mat([v[:j]for v in rows[:j]]).det()>0 for j in range(1,55))
    assert all(G[i,i]%2==0 for i in range(54))
    assert gcd(*a)==1
    label=fmpz_mat([[x]for x in a]);ar=fmpq_mat([a]);k=(ar*G.inv()*ar.transpose())[0,0]
    assert F(str(k))==F(10,3)
    H,T=label.hnf(transform=True);assert T*label==H and abs(T.det())==1
    assert H.tolist()==[[1]]+[[0]]*53
    N=fmpz_mat(T.tolist()[1:]);period=N*G*N.transpose()
    assert N*label==fmpz_mat(53,1)
    assert period.det()==10*3**24==G.det()*k
    wr=fmpz_mat([w]);assert wr*label==fmpz_mat(1,1)and(wr*G*wr.transpose())[0,0]==10
    coords=wr*T.inv();assert coords[0,0]==0
    assert all(x.denominator==1 for x in coords.tolist()[0])
    section_w=fmpz_mat([[int(x)for x in coords.tolist()[0][1:]]])
    assert section_w*N==wr and(section_w*period*section_w.transpose())[0,0]==10
    d2=F(10,4)**53/F(int(period.det()))
    assert d2==F(5**52,2**54*3**24)
    return {'status':'verified_explicit_53D_lattice','rows':[finish(53,'Explicit primitive section of Leterrier54',d2,'new_explicit_section',
            {'parent_dimension':54,'minimum_squared':10,'minimum_attained':True,'normal_label':a,
             'dual_normal_squared_norm':'10/3','period_basis':N.tolist(),'period_gram_determinant':str(period.det()),
             'minimum_witness_parent_coordinates':w,'minimum_witness_section_coordinates':section_w.tolist()[0]})],
            'source':sources['Leterrier2024'],'all_remaining_checks_exact':True,
            'parent_minimum_scope':'Published minimum10 and published computational Gram; external inputs, not independently re-enumerated'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');ap.add_argument('--source-pdf',type=Path);args=ap.parse_args();r=audit(args.source_pdf)
    if args.write_result:(P/'section53_results.json').write_text(json.dumps(r,indent=2,default=int)+'\n')
    else:assert json.loads((P/'section53_results.json').read_text())==r
    for a in r['rows']:print(a['dimension'],a['center_density'],'gain',a['gain_over_checkpoint_percent'])
