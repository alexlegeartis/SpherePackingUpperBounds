"""Exact setup audit and complete interval replay of explicit P48 extensions.

No arguments audits the exact setup and archived completed reports.
--output NEW_DIR executes all 15 full enumerations and a negative control.
--write-result with --output archives compact reports after a complete run.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
import sympy as sp
from flint import fmpz_mat, fmpq_mat

if not __debug__:
    raise RuntimeError('Run without -O.')
P = Path(__file__).resolve().parent
DATA = P/'explicit_p48'
FLAGS = ['-O3','-std=c++17','-fno-fast-math','-ffp-contract=off']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def setup():
    body = (DATA/'P48q.html').read_text().split('<a NAME="GRAM"',1)[1].split('</a><br>',1)[1].split('<a NAME=',1)[0]
    z = list(map(int,re.findall(r'-?\d+',re.sub('<[^>]*>',' ',body))))
    assert z[0] == 48 and len(z) == 1+48*49//2
    G = [[0]*48 for _ in range(48)]
    k = 1
    for i in range(48):
        for j in range(i+1):
            G[i][j] = G[j][i] = z[k]
            k += 1
    with np.load(DATA/'certificate.npz',allow_pickle=False) as d:
        assert all(np.issubdtype(d[key].dtype,np.integer) for key in d.files)
        assert d['G'].tolist() == G
        T,Q,B,B0,R = [fmpz_mat(d[key].tolist()) for key in ('T','Q','B_reduced','B_original','vertical_root_gram')]
    g = fmpz_mat(G)
    assert g.det() == 1 and abs(T.det()) == 1
    assert all(g[i,i] % 2 == 0 for i in range(48))
    assert all(fmpz_mat([v[:j] for v in G[:j]]).det() > 0 for j in range(1,49))
    assert T*g*T.transpose() == Q and Q == Q.transpose()
    assert all(int(x) % 2 == 0 for v in (B*T-B0).tolist() for x in v)
    assert R.tolist() == [[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]
    L,D = sp.Matrix(Q.tolist()).LDLdecomposition(hermitian=False)
    assert L*D*L.T == sp.Matrix(Q.tolist())
    scale, text = 2**40, []
    for i in range(48):
        assert D[i,i] > 0
        for x in [D[i,i]]+[L[j,i] for j in range(48)]:
            lo,hi = int(sp.floor(x*scale)),int(sp.ceiling(x*scale))
            assert abs(lo) < 2**53 and abs(hi) < 2**53
            assert sp.Rational(lo,scale) <= x <= sp.Rational(hi,scale)
            text.append(f'{lo} {hi}\n')
    files = {'ldl.txt':''.join(text), 'gram.txt':'\n'.join(' '.join(map(str,v)) for v in Q.tolist())+'\n'}
    labels = set()
    for bits in range(1,16):
        z = fmpz_mat([[(bits >> i) & 1 for i in range(4)]])
        c = [int(x) % 2 for x in (z*B).tolist()[0]]
        assert any(c) and tuple(c) not in labels
        labels.add(tuple(c))
        assert int((fmpz_mat([c])*Q*fmpz_mat([c]).transpose())[0,0]) % 4 == int((z*R*z.transpose())[0,0]) % 4
        files[f'coset_{bits}.txt'] = ' '.join(map(str,c))+'\n'
    # Deliberately short class: the enumerator must reject it with a vector.
    i = next(i for i in range(48) if Q[i,i] == 6)
    control = [int(j == i) for j in range(48)]
    files['negative_control.txt'] = ' '.join(map(str,control))+'\n'
    rows = []
    for r in range(1,5):
        b = fmpq_mat(B0.tolist()[:r])
        qr = fmpq_mat([v[:r] for v in R.tolist()[:r]])
        vertical = qr*5/4
        full = fmpq_mat(48+r,48+r)
        cross = fmpq_mat(g)*b.transpose()/2
        bottom = b*fmpq_mat(g)*b.transpose()/4+vertical
        for i in range(48):
            for j in range(48):full[i,j]=g[i,j]
            for j in range(r):full[i,48+j]=full[48+j,i]=cross[i,j]
        for i in range(r):
            for j in range(r):full[48+i,48+j]=bottom[i,j]
        assert full.det() == vertical.det()
        d2 = F(3,2)**(48+r)/F(str(full.det()))
        rows.append({'dimension':48+r, 'horizontal_residue_rows_original':[[int(x) for x in v] for v in B0.tolist()[:r]],
                     'vertical_gram':[[str(x) for x in v] for v in vertical.tolist()],
                     'full_gram_determinant':str(full.det()), 'minimum_squared':6,
                     'center_density_squared_exact':str(d2), 'explicit_horizontal_shifts':True})
    return files, rows, Q


def check_reports(directory, files, rows, Q):
    for name, body in files.items():
        assert (directory/name).read_text() == body, name
    manifest = json.loads((directory/'manifest.json').read_text())
    assert manifest['input_sha256'] == {name:sha(directory/name) for name in files}
    assert manifest['certificate_sha256'] == sha(DATA/'certificate.npz')
    assert manifest['catalogue_sha256'] == sha(DATA/'P48q.html')
    assert manifest['enumerator_source_sha256'] == sha(DATA/'interval_enum.cpp')
    assert manifest['compiler_flags'] == FLAGS
    reports = []
    for bits in range(1,16):
        a = json.loads((directory/f'report_{bits}.json').read_text())
        assert a['status'] == 'COMPLETE' and not a['found'] and a['norm'] == -1
        assert a['shard_id'] == 0 and a['shard_count'] == 1
        assert a['prefixes'] == a['assigned_prefixes']
        reports.append(a)
    control = json.loads((directory/'negative_control_report.json').read_text())
    assert control['status'] == 'EXACT_REJECT' and control['found']
    w = fmpz_mat([control['witness']])
    assert 0 < control['norm'] <= 10
    assert int((w*Q*w.transpose())[0,0]) == control['norm']
    expected = list(map(int, files['negative_control.txt'].split()))
    assert all((a-b) % 2 == 0 for a,b in zip(control['witness'],expected))
    return {'status':'verified_explicit_49_to_52_extensions',
            'parent':'Published P48q minimum 6; not re-enumerated',
            'positive_cosets_completed':15, 'negative_control_rejected_exactly':True,
            'total_positive_enumeration_nodes':sum(a['nodes'] for a in reports),
            'all_setup_matrix_checks_exact':True,
            'minimum_proof_arithmetic':'Outward-rounded IEEE-754 binary64 enumeration, exact dyadic input bounds, exact integer leaf checks; documented compiler/arithmetic assumptions',
            'rows':rows, 'world_record_established':False}


def run_full(out, files):
    out.mkdir(parents=True,exist_ok=False)
    for name, body in files.items():(out/name).write_text(body)
    exe = out/'interval_enum'
    subprocess.run(['g++',*FLAGS,str(DATA/'interval_enum.cpp'),'-o',str(exe)],check=True)
    manifest = {'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
                'compiler_flags':FLAGS, 'input_sha256':{name:sha(out/name) for name in files},
                'certificate_sha256':sha(DATA/'certificate.npz'), 'catalogue_sha256':sha(DATA/'P48q.html'),
                'enumerator_source_sha256':sha(DATA/'interval_enum.cpp'), 'binary_sha256':sha(exe)}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def run(label):
        negative = label == 'negative_control'
        name = label+'.txt' if negative else 'coset_'+label+'.txt'
        args = [str(exe),str(out/'ldl.txt'),str(out/'gram.txt'),str(out/name),'0','1']
        with (out/(label+'.stderr')).open('w') as log:
            p = subprocess.run(args,stdout=subprocess.PIPE,stderr=log,text=True)
        assert p.returncode == (1 if negative else 0), (label,p.stdout)
        dest = 'negative_control_report.json' if negative else 'report_'+label+'.json'
        (out/dest).write_text(json.dumps(json.loads(p.stdout),indent=2)+'\n')
        print(label,'complete',flush=True)
    run('negative_control')
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run,map(str,range(1,16))))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--output',type=Path)
    ap.add_argument('--write-result',action='store_true')
    args = ap.parse_args()
    assert not args.write_result or args.output
    files, rows, Q = setup()
    directory = args.output.resolve() if args.output else DATA/'replay'
    if args.output:run_full(directory,files)
    result = check_reports(directory,files,rows,Q)
    if args.write_result:
        dest = DATA/'replay'
        dest.mkdir(exist_ok=True)
        for name in list(files)+['manifest.json','negative_control_report.json']+[f'report_{i}.json' for i in range(1,16)]:
            shutil.copyfile(directory/name,dest/name)
        (P/'explicit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert json.loads((P/'explicit_results.json').read_text()) == result
    print(result['status'],result['positive_cosets_completed'],'cosets;',
          result['total_positive_enumeration_nodes'],'positive search nodes; negative control passed.')
    for a in result['rows']:
        print(a['dimension'],'det =',a['full_gram_determinant'],
              'minimum =',a['minimum_squared'],'explicit =',a['explicit_horizontal_shifts'])
