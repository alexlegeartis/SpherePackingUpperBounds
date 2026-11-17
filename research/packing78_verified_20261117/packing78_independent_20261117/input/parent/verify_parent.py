from pathlib import Path
import json,numpy as np
from flint import fmpz_mat,fmpq_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
D=Path(__file__).parent
def parent80():
    z=json.loads((D/'unitary_sections80.json').read_text());assert z['prime']==41 and z['field_nonsquare']==3
    field=fq_default_ctx(modulus=fmpz_mod_poly_ctx(41)([-3,0,1]),fq_type='FQ_NMOD');ring=fq_default_poly_ctx(field)
    def fe(a):
        assert isinstance(a,int)and 0<=a<41**2
        return field([a%41,a//41])
    polys=[]
    for row in z['unitaries']:
        a,b,c,d=map(fe,row);assert a*a**41+b*b**41==1 and c*c**41+d*d**41==1 and a*c**41+b*d**41==0
        x,y=-ring([b,a])**14,ring([d,c])**21;assert y*y==x*x*x+ring([1]+[0]*41+[1]);polys.append((x,y))
    count=len(polys);assert count==190;H=[[0]*count for _ in range(count)]
    for i,(x,y)in enumerate(polys):
        H[i][i]=14
        for j,(u,v)in enumerate(polys[:i]):
            dx,dy=x-u,y-v
            if dx.is_zero()and dy.is_zero():q=14
            elif dx.is_zero()and(y+v).is_zero():q=-14
            else:
                orders=[]
                if not dx.is_zero():orders.append(14-dx.degree())
                if not dy.is_zero():orders.append(21-dy.degree())
                q=7-dx.gcd(dy).degree()-min(orders)
            H[i][j]=H[j][i]=q
    with np.load(D/'parent80.npz',allow_pickle=False)as d:
        assert all(np.issubdtype(d[k].dtype,np.integer)for k in d.files)
        G,Y=[fmpz_mat(d[k].tolist())for k in ['G','section_coordinates']]
    assert G==G.transpose()and G.det()==41**12
    assert all(fmpz_mat([r[:k]for r in G.tolist()[:k]]).det()>0 for k in range(1,81))
    assert Y*G*Y.transpose()==fmpz_mat(H)
    assert Y.hnf().tolist()[:80]==[[int(i==j)for j in range(80)]for i in range(80)]
    return G,np.array(Y.tolist()[0],dtype=np.int64)


if __name__=='__main__':
 G,v=parent80();print('MW80 source verified',G.det())
