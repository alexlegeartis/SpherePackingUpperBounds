from pathlib import Path
import json,numpy as np
from flint import fmpz_mat,fmpq_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
D=Path(__file__).parent
def parent56():
    z=json.loads((D/'unitary_sections56.json').read_text());assert z['prime']==29 and z['field_nonsquare']==2
    field=fq_default_ctx(modulus=fmpz_mod_poly_ctx(29)([-2,0,1]),fq_type='FQ_NMOD');ring=fq_default_poly_ctx(field)
    def fe(a):
        assert isinstance(a,int)and 0<=a<29**2
        return field([a%29,a//29])
    polys=[]
    for row in z['unitaries']:
        a,b,c,d=map(fe,row);assert a*a**29+b*b**29==1 and c*c**29+d*d**29==1 and a*c**29+b*d**29==0
        x,y=-ring([b,a])**10,ring([d,c])**15;assert y*y==x*x*x+ring([1]+[0]*29+[1]);polys.append((x,y))
    count=len(polys);assert count==140;H=[[0]*count for _ in range(count)]
    for i,(x,y)in enumerate(polys):
        H[i][i]=10
        for j,(u,v)in enumerate(polys[:i]):
            dx,dy=x-u,y-v
            if dx.is_zero()and dy.is_zero():q=10
            elif dx.is_zero()and(y+v).is_zero():q=-10
            else:
                orders=[]
                if not dx.is_zero():orders.append(10-dx.degree())
                if not dy.is_zero():orders.append(15-dy.degree())
                q=5-dx.gcd(dy).degree()-min(orders)
            H[i][j]=H[j][i]=q
    with np.load(D/'parent56.npz',allow_pickle=False)as d:
        assert all(np.issubdtype(d[k].dtype,np.integer)for k in d.files)
        G,Y=[fmpz_mat(d[k].tolist())for k in ['G','section_coordinates']]
    assert G==G.transpose()and G.det()==29**8
    assert all(fmpz_mat([r[:k]for r in G.tolist()[:k]]).det()>0 for k in range(1,57))
    assert Y*G*Y.transpose()==fmpz_mat(H)
    assert Y.hnf().tolist()[:56]==[[int(i==j)for j in range(56)]for i in range(56)]
    return G,np.array(Y.tolist()[0],dtype=np.int64)


if __name__=='__main__':
 G,v=parent56();print('MW56 source verified',G.det())
