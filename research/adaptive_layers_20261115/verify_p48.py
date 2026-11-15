"""Exact setup and complete saved-coset audit of the new 50D-52D lattices.

--replay NEW_DIRECTORY recompiles and reruns all fifteen classes and the
negative control. Parent P48q minimum6 is an inherited published theorem.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
import argparse,json,hashlib,subprocess,sys
import numpy as np
from flint import fmpz_mat,fmpq_mat
from common import P,metric_audit,finish
if not __debug__:raise RuntimeError('Run without -O.')
D=P/'p48';SOURCE=P.parent/'covering_and_deep_cosets_20261115/explicit_p48/interval_enum.cpp'
FLAGS=['-O3','-std=c++17','-fno-fast-math','-ffp-contract=off']
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()

def setup():
    sys.path.insert(0,str(P.parent/'range48_96_20261115'))
    import verify_explicit as original
    oldfiles,_,Qold=original.setup()
    with np.load(D/'certificate.npz',allow_pickle=False)as a:
        assert all(np.issubdtype(a[k].dtype,np.integer)for k in a.files)
        G,T,Q,B,B0=[fmpz_mat(a[k].tolist())for k in('G','T','Q','B_reduced','B_original')]
    assert Q==Qold and T*G*T.transpose()==Q and abs(T.det())==G.det()==1
    assert all(int(x)%2==0 for v in(B*T-B0).tolist()for x in v)
    files={k:oldfiles[k]for k in('ldl.txt','gram.txt','negative_control.txt')}
    seen=set()
    for bits in range(1,16):
        a=[(bits>>i)&1 for i in range(4)];c=[int(x)%2 for x in(fmpz_mat([a])*B).tolist()[0]]
        assert any(c)and tuple(c)not in seen;seen.add(tuple(c))
        q=(a[2]+a[3]+a[2]*a[3])%2
        assert int((fmpz_mat([c])*Q*fmpz_mat([c]).transpose())[0,0])%4==2*q
        files[f'coset_{bits}.txt']=' '.join(map(str,c))+'\n'
    M4=[[F(2),F(-1),F(0),F(0)],[F(-1),F(2),F(-1),F(-1)],
        [F(0),F(-1),F(5,2),F(-1,4)],[F(0),F(-1),F(-1,4),F(5,2)]]
    rows=[]
    for r in range(2,5):
        M=[v[:r]for v in M4[:r]]
        targets=[F(6)if bits==0 else F(2)if bits<4 else F(5,2)for bits in range(2**r)]
        audit=metric_audit(M,targets);m=fmpq_mat([[str(x)for x in v]for v in M]);g=fmpq_mat(G)
        b=fmpq_mat(B0.tolist()[:r]);cross=g*b.transpose()/2;bottom=b*g*b.transpose()/4+m
        full=fmpq_mat(48+r,48+r)
        for i in range(48):
            for j in range(48):full[i,j]=g[i,j]
            for j in range(r):full[i,48+j]=full[48+j,i]=cross[i,j]
        for i in range(r):
            for j in range(r):full[48+i,48+j]=bottom[i,j]
        assert full.det()==m.det()
        rows.append(finish(48+r,'P48q deep plane and fitted vertical metric',F(3,2)**(48+r)/F(str(m.det())),
                           'new_explicit_lattice',{'minimum_squared':6,'gram_determinant':str(m.det()),
                           'metric_audit':audit,'explicit_shifts':True,'horizontal_residues_original':[[int(x)for x in v]for v in B0.tolist()[:r]]}))
    return files,rows,Q

def audit(directory,files,rows,Q):
    for name,text in files.items():assert(directory/name).read_text()==text
    manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['enumerator_sha256']==sha(SOURCE)and manifest['certificate_sha256']==sha(D/'certificate.npz')
    assert manifest['compiler_flags']==FLAGS
    assert manifest['input_sha256']=={name:sha(directory/name)for name in files}
    nodes=0
    for i in range(1,16):
        a=json.loads((directory/f'report_{i}.json').read_text())
        assert a['status']=='COMPLETE'and not a['found']and a['norm']==-1
        assert a['cutoff']==(12 if i<4 else 10)
        assert a['shard_id']==0 and a['shard_count']==1 and a['prefixes']==a['assigned_prefixes']
        assert a['residue']==list(map(int,files[f'coset_{i}.txt'].split()))
        nodes+=a['nodes']
    a=json.loads((directory/'negative_control_report.json').read_text())
    assert a['status']=='EXACT_REJECT'and a['found']and 0<a['norm']<=10
    w=fmpz_mat([a['witness']]);assert int((w*Q*w.transpose())[0,0])==a['norm']
    assert all((x-y)%2==0 for x,y in zip(a['witness'],map(int,files['negative_control.txt'].split())))
    return {'status':'verified_explicit_deep_plane_lattices','rows':rows,'completed_cosets':15,'completed_nodes':nodes,
            'negative_control_passed':True,'parent_minimum_input':'Published P48q minimum6; inherited',
            'enumeration_assumptions':'Outward-rounded IEEE-754 binary64; exact dyadic inputs and integer leaf checks; -fno-fast-math and -ffp-contract=off',
            'all_matrix_and_vertical_checks_exact':True}

def replay(out,files):
    out.mkdir(parents=True,exist_ok=False)
    for name,text in files.items():(out/name).write_text(text)
    for cut in(10,12):subprocess.run(['g++',*FLAGS,f'-DCUTOFF={cut}',str(SOURCE),'-o',str(out/f'enum{cut}')],check=True)
    def one(bits):
        cut=12 if 0<bits<4 else 10;name=f'coset_{bits}.txt'if bits else'negative_control.txt'
        r=subprocess.run([str(out/f'enum{cut}'),str(out/'ldl.txt'),str(out/'gram.txt'),str(out/name),'0','1'],capture_output=True,text=True)
        assert r.returncode==(0 if bits else 1)
        a=json.loads(r.stdout);a['residue']=list(map(int,files[name].split()))
        (out/(f'report_{bits}.json'if bits else'negative_control_report.json')).write_text(json.dumps(a,indent=2)+'\n')
        print('coset',bits,a['status'],flush=True)
    one(0)
    with ThreadPoolExecutor(max_workers=3)as pool:list(pool.map(one,range(1,16)))
    manifest={'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'compiler_flags':FLAGS,
              'enumerator_source':str(SOURCE.relative_to(P.parents[1])),'enumerator_sha256':sha(SOURCE),
              'certificate_sha256':sha(D/'certificate.npz'),'input_sha256':{name:sha(out/name)for name in files}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');ap.add_argument('--replay',type=Path);args=ap.parse_args()
    files,rows,Q=setup();directory=args.replay.resolve()if args.replay else D/'replay'
    if args.replay:replay(directory,files)
    result=audit(directory,files,rows,Q)
    if args.write_result:(P/'p48_results.json').write_text(json.dumps(result,indent=2)+'\n')
    else:assert json.loads((P/'p48_results.json').read_text())==result
    print(result['status'],result['completed_nodes'],'completed nodes')
    for a in rows:print(a['dimension'],a['center_density'],'gain',a['gain_over_checkpoint_percent'])
