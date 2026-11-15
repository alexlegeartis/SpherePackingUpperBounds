"""Exact MW56/CM setup and saved complete coset reports for the new 58D lattice.

--replay NEW_DIR compiles and repeats the complete positive and negative runs.
Published parent minimum10 is inherited. See PROOF.md for arithmetic assumptions.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import argparse,json,hashlib,subprocess,sys
import numpy as np
import sympy as sp
from flint import fmpz_mat,fmpq_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
if not __debug__: raise RuntimeError('Run without -O.')
P=Path(__file__).resolve().parent
FLAGS=['-O3','-std=c++17','-fno-fast-math','-ffp-contract=off']
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()

def parent():
    source=json.loads((P/'unitary_sections.json').read_text())
    assert source['prime']==29 and source['field_nonsquare']==2
    field=fq_default_ctx(modulus=fmpz_mod_poly_ctx(29)([-2,0,1]),fq_type='FQ_NMOD')
    ring=fq_default_poly_ctx(field)
    def elem(z):
        assert isinstance(z,int) and 0<=z<29**2
        return field([z%29,z//29])
    polys=[]
    for row in source['unitaries']:
        a,b,c,d=map(elem,row)
        assert a*a**29+b*b**29==1 and c*c**29+d*d**29==1 and a*c**29+b*d**29==0
        x,y=-ring([b,a])**10,ring([d,c])**15
        assert y*y==x*x*x+ring([1]+[0]*29+[1])
        polys.append((x,y))
    assert len(polys)==120
    heights=[[0]*120 for _ in range(120)]
    for i,(x,y)in enumerate(polys):
        heights[i][i]=10
        for j,(u,v)in enumerate(polys[:i]):
            dx,dy=x-u,y-v
            if dx.is_zero()and dy.is_zero():h=10
            elif dx.is_zero()and(y+v).is_zero():h=-10
            else:
                orders=[]
                if not dx.is_zero():orders.append(10-dx.degree())
                if not dy.is_zero():orders.append(15-dy.degree())
                h=5-dx.gcd(dy).degree()-min(orders)
            heights[i][j]=heights[j][i]=h
    with np.load(P/'certificate.npz',allow_pickle=False)as data:
        assert all(np.issubdtype(data[k].dtype,np.integer)for k in data.files)
        G,Y,T,Q,C=[fmpz_mat(data[k].tolist())for k in ['G','point_coordinates','T','Q','C']]
        c=[int(x)for x in data['c']]
    assert G.nrows()==G.ncols()==56 and Y.nrows()==120 and Y.ncols()==56
    assert G==G.transpose()and all(G[i,i]%2==0 for i in range(56))
    assert all(fmpz_mat([row[:k]for row in G.tolist()[:k]]).det()>0 for k in range(1,57))
    assert Y*G*Y.transpose()==fmpz_mat(heights)
    assert Y.hnf().tolist()[:56]==[[int(i==j)for j in range(56)]for i in range(56)]
    assert G.det()==29**8 and abs(T.det())==1 and T*G*T.transpose()==Q
    I=fmpz_mat(np.eye(56,dtype=int).tolist())
    assert C*Q*C.transpose()==Q and C*C+C+I==fmpz_mat(56,56)
    assert abs(C.det())==1 and any(c)and all(x in(0,1)for x in c)
    a=fmpz_mat([c]);assert int((a*Q*a.transpose())[0,0])%4==2
    return G,Q,C,a

def inputs(Q,a):
    q=sp.Matrix(Q.tolist());L,D=q.LDLdecomposition(hermitian=False)
    assert L*D*L.T==q
    lines=[];scale=2**40
    for i in range(56):
        assert D[i,i]>0
        for x in [D[i,i]]+[L[j,i]for j in range(56)]:
            lo,hi=int(sp.floor(scale*x)),int(sp.ceiling(scale*x))
            assert abs(lo)<2**53 and abs(hi)<2**53
            assert sp.Rational(lo,scale)<=x<=sp.Rational(hi,scale)
            lines.append(f'{lo} {hi}\n')
    return {'ldl.txt':''.join(lines),'gram.txt':'\n'.join(' '.join(map(str,row))for row in Q.tolist())+'\n',
            'coset.txt':' '.join(map(str,a.tolist()[0]))+'\n'}

def check_reports(directory,files,Q):
    for name,body in files.items():assert(directory/name).read_text()==body
    manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['compiler_flags']==FLAGS and manifest['cutoff']==22 and manifest['dimension']==56
    assert manifest['enumerator_sha256']==sha(P/'interval_enum.cpp')
    assert manifest['certificate_sha256']==sha(P/'certificate.npz')
    assert manifest['input_sha256']=={n:sha(directory/n)for n in [*files,'negative_control.txt']}
    counts=[]
    for name in ['report_primary.json','report_rerun.json']:
        report=json.loads((directory/name).read_text())
        assert report['status']=='COMPLETE'and not report['found']and report['norm']==-1 and report['cutoff']==22
        assert report['shard_id']==0 and report['shard_count']==1 and report['prefixes']==report['assigned_prefixes']
        counts.append(report['nodes'])
    assert counts[0]==counts[1]==2904500151
    negative=json.loads((directory/'negative_control_report.json').read_text())
    assert negative['status']=='EXACT_REJECT'and negative['found']and 0<negative['norm']<=22 and negative['cutoff']==22
    v=fmpz_mat([negative['witness']]);assert int((v*Q*v.transpose())[0,0])==negative['norm']
    residue=list(map(int,(directory/'negative_control.txt').read_text().split()))
    assert all((int(x)-y)%2==0 for x,y in zip(v.tolist()[0],residue))
    return {'completed_nodes_per_run':counts[0],'complete_full_runs':2,'negative_control_passed':True}

def audit():
    G,Q,C,a=parent();files=inputs(Q,a);enumeration=check_reports(P/'replay',files,Q)
    b=fmpz_mat([[int(x)%2 for x in(a*C).tolist()[0]]]);ab=fmpz_mat([[int(x)%2 for x in(a+b).tolist()[0]]])
    assert all((int(x)-int(y))%2==0 for x,y in zip((a*C*C).tolist()[0],ab.tolist()[0]))
    parities={tuple(a.tolist()[0]),tuple(b.tolist()[0]),tuple(ab.tolist()[0])}
    assert len(parities)==3 and all(any(row)for row in parities)
    for row in parities:
        v=fmpz_mat([list(row)]);assert int((v*Q*v.transpose())[0,0])%4==2
    B=fmpq_mat([a.tolist()[0],b.tolist()[0]])
    M=fmpq_mat([['7/2','-7/4'],['-7/4','7/2']]);assert str(M.det())=='147/16'
    # z1²-z1*z2+z2² is a positive integer for each nonzero integer pair.
    # Thus M has minimum7/2; all-even nonzero coefficients have norm>=14.
    full=fmpq_mat(58,58);cross=fmpq_mat(Q)*B.transpose()/2
    bottom=B*fmpq_mat(Q)*B.transpose()/4+M
    for i in range(56):
        for j in range(56):full[i,j]=Q[i,j]
        for j in range(2):full[i,56+j]=full[56+j,i]=cross[i,j]
    for i in range(2):
        for j in range(2):full[56+i,56+j]=bottom[i,j]
    assert full.det()==G.det()*M.det()
    witness=next(i for i in range(56)if Q[i,i]==10)
    d2=F(5,2)**58/F(str(full.det()))
    assert d2==F(34694469519536141888238489627838134765625,1324710820694697074011278409728)
    comparison=F(2**36,3)
    external=json.loads((P.parent/'novel_routes_20261113/audits/gamma72_leech_sections_verified.json').read_text())
    old=next(x for x in external['sections']if x['dimension']==58)
    assert external['all_passed']and F(old['center_density_squared'])==comparison
    assert int(old['normal_determinant'])==3*2**22
    ratio2=d2/comparison
    ratio=F(10*5**28,7*2**46*29**4);assert ratio*ratio==ratio2
    with localcontext()as ctx:
        ctx.prec=75;density=(Decimal(d2.numerator)/d2.denominator).sqrt();gain=100*((Decimal(ratio2.numerator)/ratio2.denominator).sqrt()-1);k=int(gain*10**6)
    lo,hi=F(k,10**6),F(k+1,10**6)
    assert(1+lo/100)**2<=ratio2<(1+hi/100)**2
    # Secondary directly reported comparator: Vardy1999 p.278 gives
    # Elkies's58D log2density prefix16.46. Use its upper endpoint16.47.
    with localcontext()as ctx:
        ctx.prec=75
        reported_gain=100*((Decimal(d2.numerator)/d2.denominator).sqrt()/Decimal(2)**Decimal('16.47')-1)
        reported_floor=F(int(reported_gain*10**6),10**6)
        reported_display=f'>{Decimal(reported_floor.numerator)/reported_floor.denominator:.6f}%'
    assert d2**50>2**1647*(1+reported_floor/100)**100
    result={'status':'verified_explicit_58D_CM_lattice','all_setup_and_layer_checks_exact':True,
        'maximum_parent_or_target_dimension':72,'parent':{'dimension':56,'minimum_squared_inherited':10,'gram_determinant':str(G.det()),'polynomial_sections_checked':120,'direct_polynomial_gcd_heights_checked':7140,'full_integer_span_checked':True},
        'secondary_reported_comparison':{'source':'Elkies guarantee reported by Vardy1999, p.278','log2_density_prefix':'16.46','prefix_upper_endpoint':'16.47','gain_display':reported_display,'exact_integer_power_comparison_passed':True},
        'enumeration':enumeration,'arithmetic_assumptions':'Outward-rounded IEEE-754 binary64 enumeration; exact dyadic LDL bounds and integer leaf checks; -fno-fast-math -ffp-contract=off. Matrix, field and layer checks are exact.',
        'rows':[{'dimension':58,'name':'Explicit MW56 CM plane with one certified deep coset','kind':'new_explicit_lattice','density_squared':str(d2),'center_density':str(density),'minimum_squared':10,'minimum_attained_in_parent_basis_index':witness,'gram_determinant':str(full.det()),'horizontal_residues':B.tolist(),'CM_transports_all_three_nonzero_classes':True,'half_coset_squared_depth_at_least':'13/2','vertical_gram':M.tolist(),'comparison_density_squared':str(comparison),'comparison_name':'Classical Gamma72/Leech-block rank14 section, derived from Nebe2012','ratio_exact':str(ratio),'ratio_squared':str(ratio2),'gain_percent':str(gain),'gain_interval':[str(lo),str(hi)],'gain_display':f'{gain:.6f}%'}]}
    return result,full

def replay(out):
    out.mkdir(parents=True,exist_ok=False)
    for name in ['ldl.txt','gram.txt','coset.txt','negative_control.txt']:(out/name).write_bytes((P/'replay'/name).read_bytes())
    exe=out/'enum22';subprocess.run(['g++',*FLAGS,'-DCUTOFF=22',str(P/'interval_enum.cpp'),'-o',str(exe)],check=True)
    for label,source,code in [('negative_control_report.json','negative_control.txt',1),('report_primary.json','coset.txt',0),('report_rerun.json','coset.txt',0)]:
        with(out/(label+'.stderr')).open('w')as err:
            run=subprocess.run([str(exe),str(out/'ldl.txt'),str(out/'gram.txt'),str(out/source),'0','1'],stdout=subprocess.PIPE,stderr=err,text=True)
        assert run.returncode==code
        (out/label).write_text(run.stdout+'\n');print(label,'completed',flush=True)
    (out/'manifest.json').write_bytes((P/'replay/manifest.json').read_bytes())
    _,Q,_,a=parent();return check_reports(out,inputs(Q,a),Q)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-result',action='store_true');parser.add_argument('--replay',type=Path);args=parser.parse_args()
    result,gram=audit()
    if args.replay:assert replay(args.replay)==result['enumeration']
    if args.write_result:
        (P/'results.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
        (P/'gram58.json').write_text(json.dumps([[str(x)for x in row]for row in gram.tolist()],indent=2)+'\n')
    else:
        assert json.loads((P/'results.json').read_text())==json.loads(json.dumps(result,default=str))
        assert json.loads((P/'gram58.json').read_text())==[[str(x)for x in row]for row in gram.tolist()]
    row=result['rows'][0];print('58D',row['center_density'],'gain',row['gain_display'],'two completed runs',result['enumeration']['completed_nodes_per_run'])
