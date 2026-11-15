"""Exact local geometry and the counterexample rejecting the 92D deformation.

Only58D is a new certified density improvement. Local optimality below means
local optimality among lattice Gram deformations, not global optimality.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json,argparse,importlib.util
import numpy as np
from flint import fmpz_mat,fmpq_mat,nmod_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
from verify import parent as parent56
if not __debug__:raise RuntimeError('Run without -O.')
P=Path(__file__).resolve().parent;D=P/'frontier'

def exact_norms(V,G,target):
    assert len(V)>0 and np.issubdtype(V.dtype,np.integer)and np.issubdtype(G.dtype,np.integer)
    bound=G.shape[0]**2*int(abs(V).max())**2*int(abs(G).max())
    assert bound<2**63
    assert np.all(np.einsum('vi,ij,vj->v',V,G,V)==target)

def dyad_rank(V):
    i,j=np.triu_indices(V.shape[1])
    return nmod_mat((V[:,i]*V[:,j]).tolist(),65521).rank()

def orbit(G,seed,maps,target):
    assert np.issubdtype(maps.dtype,np.integer)
    n=G.nrows();gg=np.array(G.tolist(),dtype=np.int64)
    for C in maps:
        c=fmpz_mat(C.tolist());assert c*G*c.transpose()==G and abs(c.det())==1
    def key(v):return tuple(map(int,v if v[np.flatnonzero(v)[0]]>0 else-v))
    rows=[np.array(key(seed),dtype=np.int64)];seen={tuple(rows[0])};pos=0
    while pos<len(rows):
        v=rows[pos];pos+=1
        for C in maps:
            w=v@C;k=key(w)
            if k not in seen:seen.add(k);rows.append(np.array(k,dtype=np.int64))
        assert len(rows)<=10000
    V=np.array(rows);exact_norms(V,gg,target)
    moment=fmpz_mat(V.tolist()).transpose()*fmpz_mat(V.tolist())*G
    scalar=moment[0,0];assert scalar>0
    assert moment==scalar*fmpz_mat(np.eye(n,dtype=int).tolist())
    return V,str(scalar)

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
    count=len(polys);assert count==160;H=[[0]*count for _ in range(count)]
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

def parents54_55():
    raw=json.loads((D/'polynomial_points54.json').read_text());assert raw['field_modulus']==[2,2,1,0,2,0,1]
    field=fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)(raw['field_modulus']),fq_type='FQ_NMOD');ring=fq_default_poly_ctx(field);t=ring([0,1])
    def fe(a):
        assert isinstance(a,int)and 0<=a<729
        digits=[]
        for _ in range(6):digits.append(a%3);a//=3
        return field(digits)
    polys=[]
    for a in raw['points']:
        x,y=ring(list(map(fe,a['x']))),ring(list(map(fe,a['y'])))
        assert y*y==x*x*x+x+t**28 and x.degree()==10
        polys.append((x,y))
    def canonical(num,den):return F(max(9*num.degree(),84+9*den.degree()),9)
    def h(x,y,u,v,diag1,diag2):
        if x==u:
            if y==v:return diag1
            assert y==-v;return-diag1
        den=(u-x)**2;num=(v-y)**2-(x+u)*den;g=num.gcd(den);num,den=num//g,den//g
        return(canonical(num,den)-diag1-diag2)/2
    count=len(polys);H=[[0]*count for _ in range(count)]
    for i,(x,y)in enumerate(polys):
        H[i][i]=10
        for j,(u,v)in enumerate(polys[:i]):
            q=h(x,y,u,v,F(10),F(10));assert q.denominator==1;H[i][j]=H[j][i]=q.numerator
    with np.load(D/'local54_55.npz',allow_pickle=False)as d:
        assert all(np.issubdtype(d[k].dtype,np.integer)for k in d.files)
        G,Y,q3,roots,closest3,G55=[fmpz_mat(d[k].tolist())for k in ['G54','Y','q3','roots54','closest3','G55']]
    assert G==G.transpose()and G.det()==3**25
    assert all(fmpz_mat([r[:k]for r in G.tolist()[:k]]).det()>0 for k in range(1,55))
    assert Y*G*Y.transpose()==fmpz_mat(H)
    assert Y.hnf().tolist()[:54]==[[int(i==j)for j in range(54)]for i in range(54)]
    q=fmpq_mat(q3)/3
    hQ=fmpq_mat([[str(h(x,y,ring(0),t**14,F(10),F(28,3)))]for x,y in polys])
    assert fmpq_mat(Y*G)*q.transpose()==hQ
    assert str((q*fmpq_mat(G)*q.transpose())[0,0])=='28/3'
    assert any(x.denominator==3 for row in q.tolist()for x in row)
    expected=fmpq_mat(55,55);label=fmpq_mat(G)*q.transpose()
    for i in range(54):
        for j in range(54):expected[i,j]=G[i,j]
        expected[i,54]=expected[54,i]=label[i,0]
    expected[54,54]=10
    assert expected==fmpq_mat(G55)and G55.det()==2*3**24
    V=np.array(roots.tolist(),dtype=np.int64);exact_norms(V,np.array(G.tolist(),dtype=np.int64),10)
    assert roots.transpose()*roots*G==1620*fmpz_mat(np.eye(54,dtype=int).tolist())
    r54=dyad_rank(V);assert r54==1485
    W=fmpq_mat(closest3)/3;Wn=np.array(closest3.tolist(),dtype=np.int64)
    exact_norms(Wn,np.array(G.tolist(),dtype=np.int64),84)
    assert W.rank()==54 and all(x==0 for row in(fmpq_mat([[1]*W.nrows()])*W).tolist()for x in row)
    assert W.transpose()*W*fmpq_mat(G)==378*fmpq_mat(np.eye(54,dtype=int).tolist())
    layers=[]
    for row in W.tolist():
        v=fmpq_mat([row])-q;assert all(x.denominator==1 for r in v.tolist()for x in r)
        layers.append([int(x)for x in v.tolist()[0]]+[1])
    root55=fmpq_mat([row+[0]for row in roots.tolist()]);L=fmpq_mat(layers)
    assert root55.transpose()*root55/2187+L.transpose()*L/1458==fmpq_mat(G55).inv()
    full=np.array([list(map(int,row))+[0]for row in roots.tolist()]+layers,dtype=np.int64)
    exact_norms(full,np.array(G55.tolist(),dtype=np.int64),10);r55=dyad_rank(full);assert r55==1540
    return {'dimensions':[54,55],'status':'certified_local_optimality_among_lattice_Gram_deformations',
        'minimum_input':'Leterrier2024 Theorem3.2.7, Proposition3.2.22 and Theorem3.4.1; minimum10',
        'source_polynomial_points':count,'minimal_vectors54':roots.nrows(),'closest_vectors_for55':W.nrows(),
        'perfection_ranks':[r54,r55],'positive_eutaxy_weights55':['1/2187','1/1458'],
        'global_optimality_claimed':False}

def audit92():
    path=P.parent/'mw92_sections_20261115/verify.py';spec=importlib.util.spec_from_file_location('mw92_source_audit',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);G,report=m.parent()
    with np.load(P.parent/'mw92_sections_20261115/parent_certificate.npz',allow_pickle=False)as d:seed=d['section_coordinates'][0]
    maps=np.load(D/'isometries92.npz',allow_pickle=False)['maps'];V,scalar=orbit(G,seed,maps,16);rank=dyad_rank(V);assert rank==4276
    forms=np.load(D/'null_forms92.npz',allow_pickle=False)['H'];gg=np.array(G.tolist(),dtype=np.int64);I=fmpq_mat(np.eye(92,dtype=int).tolist())
    assert np.issubdtype(forms.dtype,np.integer)
    for raw in forms:
        H=fmpz_mat(raw.tolist());assert H==H.transpose();A=G.inv()*H
        assert all(x.denominator==1 for row in A.tolist()for x in row)
        assert A*A==3*I and sum(A[i,i]for i in range(92))==0
        assert 92**2*int(abs(V).max())**2*int(abs(raw).max())<2**63
        assert np.all(np.einsum('vi,ij,vj->v',V,raw,V)==0)
    h0,h1=[fmpz_mat(a.tolist())for a in forms]
    assert any(h0[i,j]*h1[k,l]!=h0[k,l]*h1[i,j]for i,j,k,l in [(0,1,0,2)])
    candidate=fmpq_mat(G)+fmpq_mat(h0)/16
    assert F(str(candidate.det()))==F(str(G.det()))*(F(253,256)**46)
    # Spectral bound: A²=3I, traceA=0. Since3/256<1/81,
    # G+H/16>(8/9)G. All parent vectors of norm>=18 stay above16.
    assert F(3,256)<F(1,81)
    # The missing shortest vectors are now explicit. CM rotates their two
    # null-form values into a triangle containing zero in its interior.
    raw=json.loads((D/'counterexample92.json').read_text());w=fmpz_mat([raw['witness']])
    assert w.nrows()==1 and w.ncols()==92 and (w*G*w.transpose())[0,0]==16
    A0=h0*G.inv();A1=h1*G.inv();C=A0*A1/3
    assert all(x.denominator==1 for row in C.tolist()for x in row)
    assert C*fmpq_mat(G)*C.transpose()==fmpq_mat(G) and C*C+C+I==0*I
    witnesses=[fmpq_mat(w),fmpq_mat(w)*C,fmpq_mat(w)*C*C]
    values=[]
    for v in witnesses:
        assert all(x.denominator==1 for row in v.tolist()for x in row)
        assert (v*fmpq_mat(G)*v.transpose())[0,0]==16
        values.append([int((v*fmpq_mat(h)*v.transpose())[0,0])for h in (h0,h1)])
    assert values[0]==[-24,24]
    assert all(sum(v[j]for v in values)==0 for j in range(2))
    assert fmpz_mat(values).rank()==2
    assert str((fmpq_mat(w)*candidate*fmpq_mat(w).transpose())[0,0])=='29/2'
    # The original positive-eutaxy shell forces trace(G^-1 D)>=0 for every
    # admissible deformation. Equality makes D a combination of H0,H1.
    # The three new vectors' positive-spanning values then force D=0.
    # This proves strict local determinant minimality; see FRONTIER.md.
    return {'dimension':92,'status':'REJECTED_BY_EXACT_SHORTEST_VECTOR',
        'known_minimal_rays':len(V),'exact_known_dyad_rank':rank,'missing_form_directions':2,
        'rejected_candidate':'G+H0/16','determinant_ratio':str(F(253,256)**46),
        'all_parent_vectors_of_norm_at_least18_proved_safe':True,
        'counterexample_parent_norm':16,'counterexample_H0_value':-24,
        'counterexample_deformed_norm':'29/2',
        'three_new_minimal_rays_null_values':values,'extended_perfection_rank':4278,
        'local_optimality_among_lattice_Gram_deformations':True,
        'registered_as_SOTA':False}

def audit():
    out=[];out.append(parents54_55());print('54D/55D local geometry verified',flush=True)
    G,_,_,_=parent56();seed=np.load(P/'certificate.npz',allow_pickle=False)['point_coordinates'][0]
    V,scalar=orbit(G,seed,np.load(D/'isometries56.npz',allow_pickle=False)['maps'],10);rank=dyad_rank(V);assert rank==1596
    out.append({'dimension':56,'minimal_rays':len(V),'perfection_rank':rank,'eutaxy_scalar':scalar,'status':'certified_local_optimality_among_lattice_Gram_deformations'});print('56D local geometry verified',flush=True)
    G,seed=parent80();V,scalar=orbit(G,seed,np.load(D/'isometries80.npz',allow_pickle=False)['maps'],14);rank=dyad_rank(V);assert rank==3240
    out.append({'dimension':80,'minimal_rays':len(V),'perfection_rank':rank,'eutaxy_scalar':scalar,'status':'certified_local_optimality_among_lattice_Gram_deformations'});print('80D local geometry verified',flush=True)
    for n,minimum in [(48,6),(72,8)]:
        assert 12*(minimum//2)>n//2+4
        out.append({'dimension':n,'status':'classical_local_optimality_from_harmonic_theta_vanishing',
            'argument':'Degree2 and4 harmonic theta series vanish by their initial order and the level1 valence bound; the minimum shell is a4-design.',
            'minimum_inherited':minimum,'global_optimality_claimed':False})
    ceiling=F(3**32,4)*F(63,64)**63;assert ceiling**25<2**1196
    out.append({'dimension':63,'status':'Ne64_inherited_minimum_antipode_bound_ruled_out',
        'universal_density_squared_ceiling':str(ceiling),
        'comparison':'Elkies guarantee reported by Vardy1999, log2center density starts23.92',
        'scope':'The inherited-minimum12 guarantee for every primitive normal and every cluster with at least2 integer labels. A stronger actual minimum would require a separate proof.'})
    # Exhaust every finite cluster in the scaled-E8 projection belonging
    # to a4E8 normal, using its half-integer norm spectrum and residue bounds.
    certificate=P.parent/'novel_routes_20261113/higher_parent/gamma72_root_cluster_rank8.npz'
    K=fmpz_mat(np.load(certificate,allow_pickle=False)['K'].tolist())
    E=fmpq_mat(K)/4
    assert E.det()==1 and all(x.denominator==1 for row in E.tolist()for x in row)
    assert all(E[i,i]==2 for i in range(8))
    assert all(fmpq_mat([[str(E[i,j])for j in range(k)]for i in range(k)]).det()>0 for k in range(1,9))
    assert all(x.denominator==1 for row in E.inv().tolist()for x in row)
    inverse=E.inv();assert all(int(inverse[i,i])%2==0 for i in range(8))
    # beta=1/2: equilateral simplex, M<=9.
    # beta=1: two-distance polynomials are independent in a54D function space.
    # Remaining beta ranges: injection modulo2,3,4 bounds M by2^8,3^8,4^8.
    factors=[F(1),9*F(15,16)**32,54*F(7,8)**32,
             2**8*F(13,16)**32,3**8*F(3,4)**32,4**8*F(7,16)**32]
    maximum=max(factors);assert maximum==9*F(15,16)**32
    density_ceiling=2**24*maximum;assert density_ceiling<3**16
    out.append({'dimension':64,'status':'Gamma72_fourE8_inherited_minimum_antipode_family_ruled_out',
        'all_cluster_sizes_and_diameters_covered':True,
        'density_guarantee_ceiling':str(density_ceiling),
        'piecewise_ratio_upper_bounds':[str(x)for x in factors],
        'scope':'Fixed4E8 normal, every finite integer-label cluster, inherited minimum8 only. Other normals or a stronger actual minimum remain possible.'})
    rejected=audit92();print('92D deformation rejected by exact witness; local optimality proved',flush=True)
    return {'review_date':'2026-11-16','new_certified_density_dimensions':[58],
        'local_and_route_audits':out,'rejected_candidates':[rejected],
        'maximum_parent_or_target_dimension':92,'new_SOTA_claims_from_frontier_audit':False}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-result',action='store_true');args=parser.parse_args();result=audit()
    if args.write_result:(P/'frontier_results.json').write_text(json.dumps(result,indent=2)+'\n')
    else:assert json.loads((P/'frontier_results.json').read_text())==result
    print('All frontier checks passed; only58D is registered as a new density improvement.')
