"""Full replay and exact setup audit for improved explicit 49--52D lattices.

One coset is enumerated through norm 12; the other fourteen through 10.
No arguments audit archived completed reports. --output NEW_DIR reruns all.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,argparse,subprocess,sys,shutil
import numpy as np
from flint import fmpz_mat,fmpq_mat

if not __debug__:
    raise RuntimeError('Run without -O.')
P=Path(__file__).resolve().parent;DATA=P/'explicit_p48'
sys.path.insert(0,str(P.parent/'range48_96_20261115'))
import verify_explicit as old
sys.path.insert(0,str(P.parent/'cross_dimension_lower_20261115'))
from verify_counting import decsqrt
FLAGS=['-O3','-std=c++17','-fno-fast-math','-ffp-contract=off']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def setup():
    original,_,Qold=old.setup()
    with np.load(DATA/'certificate.npz',allow_pickle=False) as d:
        assert all(np.issubdtype(d[k].dtype,np.integer) for k in d.files)
        G,T,Q,B,B0=[fmpz_mat(d[k].tolist()) for k in ('G','T','Q','B_reduced','B_original')]
    assert Q==Qold and T*G*T.transpose()==Q
    assert abs(T.det())==1 and G.det()==1
    assert all(int(x)%2==0 for v in (B*T-B0).tolist() for x in v)
    files={k:original[k] for k in ('ldl.txt','gram.txt','negative_control.txt')}
    seen=set()
    for bits in range(1,16):
        a=[(bits>>i)&1 for i in range(4)];z=fmpz_mat([a])
        c=[int(x)%2 for x in (z*B).tolist()[0]]
        assert any(c) and tuple(c) not in seen;seen.add(tuple(c))
        q=(a[1]+a[2]+a[1]*a[2])%2
        assert int((fmpz_mat([c])*Q*fmpz_mat([c]).transpose())[0,0])%4==2*q
        files[f'coset_{bits}.txt']=' '.join(map(str,c))+'\n'
    M4=[[F(2),F(1),F(1),F(0)],
        [F(1),F(5,2),F(5,4),F(-3,2)],
        [F(1),F(5,4),F(5,2),F(-3,2)],
        [F(0),F(-3,2),F(-3,2),F(3)]]
    rows=[]
    for r in range(1,5):
        M=[v[:r] for v in M4[:r]]
        m=fmpq_mat([[str(x) for x in v] for v in M])
        assert all(fmpq_mat([[str(x)for x in v[:j]]for v in M[:j]]).det()>0 for j in range(1,r+1))
        inv=m.inv();maximum=max(F(str(inv[i,i]))*6 for i in range(r));bound=1
        while bound*bound<=maximum:bound+=1
        # If zMz^T < 6, each |z_i|^2 < 6*(M^-1)_ii. Thus this
        # finite box contains every vector which could violate a layer bound.
        assert bound*bound>maximum
        gi=np.array([[int(4*x) for x in v]for v in M],dtype=np.int64)
        for v in product(range(-bound,bound+1),repeat=r):
            if not any(v):continue
            a=list(v)+[0]*(4-r);bits=sum((x%2)<<i for i,x in enumerate(v))
            q=(a[1]+a[2]+a[1]*a[2])%2
            target=F(6) if bits==0 else F(2) if bits==1 else F(5,2) if q else F(3)
            z=np.array(v,dtype=np.int64);norm=F(int(z@gi@z),4)
            assert norm>=target,(r,v,norm,target)
        b=fmpq_mat(B0.tolist()[:r]);g=fmpq_mat(G)
        full=fmpq_mat(48+r,48+r);cross=g*b.transpose()/2;bottom=b*g*b.transpose()/4+m
        for i in range(48):
            for j in range(48):full[i,j]=g[i,j]
            for j in range(r):full[i,48+j]=full[48+j,i]=cross[i,j]
        for i in range(r):
            for j in range(r):full[48+i,48+j]=bottom[i,j]
        assert full.det()==m.det()
        d2=F(3,2)**(48+r)/F(str(m.det()))
        previous=json.loads((P/'previous_bounds.json').read_text())
        before=F(next(a for a in previous['rows']if a['dimension']==48+r)['density_squared'])
        assert d2>before
        rows.append({'dimension':48+r,'horizontal_residues_original':[[int(x)for x in v]for v in B0.tolist()[:r]],
                     'vertical_gram':[[str(x)for x in v]for v in M],'gram_determinant':str(m.det()),
                     'minimum_squared':6,'density_squared':str(d2),'center_density':decsqrt(d2),
                     'explicit_shifts':True,'small_vertical_enumeration_box':bound,
                     'ratio_to_checkpoint_squared_exact':str(d2/before)})
    return files,rows,Q


def check(directory,files,rows,Q):
    for name,body in files.items():assert (directory/name).read_text()==body
    manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['input_sha256']=={name:sha(directory/name) for name in files}
    assert manifest['certificate_sha256']==sha(DATA/'certificate.npz')
    assert manifest['enumerator_source_sha256']==sha(DATA/'interval_enum.cpp')
    assert manifest['compiler_flags']==FLAGS
    nodes=0
    for bits in range(1,16):
        a=json.loads((directory/f'report_{bits}.json').read_text())
        assert a['status']=='COMPLETE' and not a['found'] and a['norm']==-1
        assert a['cutoff']==(12 if bits==1 else 10)
        assert a['shard_id']==0 and a['shard_count']==1 and a['prefixes']==a['assigned_prefixes']
        nodes+=a['nodes']
    a=json.loads((directory/'negative_control_report.json').read_text())
    assert a['status']=='EXACT_REJECT' and a['found'] and 0<a['norm']<=10
    w=fmpz_mat([a['witness']]);assert int((w*Q*w.transpose())[0,0])==a['norm']
    assert all((x-y)%2==0 for x,y in zip(a['witness'],map(int,files['negative_control.txt'].split())))
    return {'status':'verified_deeper_explicit_49_to_52','rows':rows,'positive_cosets_completed':15,
            'positive_nodes':nodes,'deep_coset_squared_depth_at_least':'4',
            'other_singular_squared_depth_at_least':'3','nonsingular_squared_depth_at_least':'7/2',
            'negative_control_passed':True,'all_matrix_and_vertical_checks_exact':True,
            'coset_enumeration_arithmetic':'Outward-rounded IEEE-754 binary64, exact dyadic inputs and integer leaf checks; same documented assumptions as the preceding P48 verifier',
            'parent_minimum_input':'Published P48q minimum 6; not re-enumerated','worldwide_priority_proved':False}


def run(out,files):
    out.mkdir(parents=True,exist_ok=False)
    for name,body in files.items():(out/name).write_text(body)
    for cutoff in (10,12):
        subprocess.run(['g++',*FLAGS,f'-DCUTOFF={cutoff}',str(DATA/'interval_enum.cpp'),'-o',str(out/f'enum{cutoff}')],check=True)
    meta={'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
          'compiler_flags':FLAGS,'cutoffs':[10,12],
          'certificate_sha256':sha(DATA/'certificate.npz'),'enumerator_source_sha256':sha(DATA/'interval_enum.cpp'),
          'input_sha256':{name:sha(out/name)for name in files},
          'binary_sha256':{str(c):sha(out/f'enum{c}')for c in (10,12)}}
    (out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    def one(bits):
        cutoff=12 if bits==1 else 10
        name=f'coset_{bits}.txt' if bits else 'negative_control.txt'
        p=subprocess.run([str(out/f'enum{cutoff}'),str(out/'ldl.txt'),str(out/'gram.txt'),str(out/name),'0','1'],capture_output=True,text=True)
        assert p.returncode==(0 if bits else 1),(bits,p.stdout,p.stderr)
        dest=f'report_{bits}.json'if bits else 'negative_control_report.json'
        (out/dest).write_text(json.dumps(json.loads(p.stdout),indent=2)+'\n');print(bits,'complete',flush=True)
    one(0)
    with ThreadPoolExecutor(max_workers=2)as pool:list(pool.map(one,range(1,16)))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);ap.add_argument('--write-result',action='store_true');args=ap.parse_args()
    assert not args.write_result or args.output
    files,rows,Q=setup();directory=args.output.resolve()if args.output else DATA/'replay'
    if args.output:run(directory,files)
    result=check(directory,files,rows,Q)
    if args.write_result:
        dest=DATA/'replay';dest.mkdir(exist_ok=True)
        for name in list(files)+['manifest.json','negative_control_report.json']+[f'report_{i}.json'for i in range(1,16)]:shutil.copyfile(directory/name,dest/name)
        (P/'deeper_results.json').write_text(json.dumps(result,indent=2)+'\n')
    else:assert json.loads((P/'deeper_results.json').read_text())==result
    print(result['status'],result['positive_nodes'],'nodes')
    for a in rows:print(a['dimension'],a['center_density'],a['gram_determinant'])
