"""Arb certificates for a general spherical polynomial below the V7 far kernel.

The subdivision tree format and root tiling follow Felpix-Studios/sphere-packing
(Apache-2.0), commit bb97550ba47ad92bacfee567843888be16895a23,
preprints/D4/verify/stage5_kernelbox/src/kbox_run.py and check_kernelbox.py.
Modified 15 November 2026: arbitrary rational Gegenbauer polynomials, explicit
failure exceptions, one shared producer/replay interface, and bounded search.
The unchanged upstream interval engine and kernel remain in vendor/D4/.
"""
from pathlib import Path
from fractions import Fraction as Q
import sys, json, argparse, multiprocessing as mp, time

HERE=Path(__file__).resolve().parent
SRC=HERE/'vendor/D4/verify/stage5_kernelbox/src'
sys.path.insert(0,str(SRC))
import kcert as kc
import kformula as kf
import kbox_run as kr

A=Q(133,250)
MIN_WIDTH=Q(1,10**9)


def polynomial_box(t,coefficients):
    u0=kc.T.const(kf.A(1))
    value=u0.scale(kf.A(coefficients[0]))
    if len(coefficients)==1: return value
    u1=2*t
    if coefficients[1]: value=value+u1.scale(kf.A(coefficients[1]/2))
    for d in range(2,len(coefficients)):
        u0,u1=u1,2*t*u1-u0
        if coefficients[d]: value=value+u1.scale(kf.A(coefficients[d]/(d+1)))
    return value


def difference(box,coefficients,open_r,open_s):
    r,s,t=kc.box_vars(box)
    lower=kc.minusK_lb(r,s,t,kc.dec(box[0],open_r),kc.dec(box[1],open_s))
    return lower-polynomial_box(t,coefficients)


def produce_one(task):
    root,coefficients=task
    br,bs,bt,open_r,open_s=root
    count=0
    leaves=0
    minimum=None
    def visit(box,depth):
        nonlocal count,leaves,minimum
        count+=1
        if count>200000 or depth>90: raise ValueError('subdivision budget exceeded')
        f=difference(box,coefficients,open_r,open_s)
        lo=f.lower()
        if lo>=0:
            leaves+=1
            value=float(lo.mid())
            minimum=value if minimum is None else min(minimum,value)
            return 'L'
        candidates=[d for d in range(3) if box[d][1]-box[d][0]>2*MIN_WIDTH]
        if not candidates: raise ValueError('unresolved box: '+str(box))
        scores=kc.split_scores(f)
        d=max(candidates,key=lambda i:scores[i])
        left=list(box);right=list(box)
        lo,hi=box[d];mid=(lo+hi)/2
        left[d]=(lo,mid);right[d]=(mid,hi)
        return [d,visit(tuple(left),depth+1),visit(tuple(right),depth+1)]
    tree=visit((br,bs,bt),0)
    return dict(r=list(map(str,br)),s=list(map(str,bs)),t=list(map(str,bt)),
                open_r=open_r,open_s=open_s,tree=tree),dict(evaluations=count,leaves=leaves,minimum=minimum)


def check_tree(path,j,coefficients,workers):
    data=json.loads(Path(path).read_text())
    if data['j']!=j or list(map(Q,data['coefficients']))!=coefficients or Q(data['a'])!=A:
        raise ValueError('certificate parameters differ')
    if type(data['nt']) is not int or not 1<=data['nt']<=10000:
        raise ValueError('invalid root chunk count')
    expected=kr.roots(j,data['nt'])
    encoded=data['roots']
    got=[(tuple(map(Q,e['r'])),tuple(map(Q,e['s'])),tuple(map(Q,e['t'])),e['open_r'],e['open_s']) for e in encoded]
    if any(type(e['open_r']) is not bool or type(e['open_s']) is not bool for e in encoded):
        raise ValueError('radius endpoint flags must be booleans')
    if got!=expected: raise ValueError('roots do not tile the full domain')
    with mp.Pool(workers) as pool:
        # Explicitly pass trees separately; neither reported minima nor producer
        # statuses are trusted during the interval replay.
        result=pool.map(replay_entry,[(e,coefficients) for e in encoded],chunksize=4)
    return dict(roots=len(encoded),leaves=sum(x['leaves'] for x in result),
                minimum_leaf_lower=min(x['minimum'] for x in result))


def replay_entry(task):
    encoded,coefficients=task
    root=(tuple(map(Q,encoded['r'])),tuple(map(Q,encoded['s'])),tuple(map(Q,encoded['t'])),
          encoded['open_r'],encoded['open_s'])
    leaves=0
    minimum=None
    def walk(box,node,depth):
        nonlocal leaves,minimum
        if depth>100: raise ValueError('tree too deep')
        if node=='L':
            lower=difference(box,coefficients,root[3],root[4]).lower()
            if not lower>=0: raise ValueError('uncertified leaf: '+str(box))
            leaves+=1
            minimum=float(lower.mid()) if minimum is None else min(minimum,float(lower.mid()))
            return
        if not (type(node) is list and len(node)==3 and type(node[0]) is int and node[0] in (0,1,2)):
            raise ValueError('invalid tree node')
        d=node[0];lo,hi=box[d];mid=(lo+hi)/2
        left=list(box);right=list(box)
        left[d]=(lo,mid);right[d]=(mid,hi)
        walk(tuple(left),node[1],depth+1)
        walk(tuple(right),node[2],depth+1)
    walk(root[:3],encoded['tree'],0)
    return dict(leaves=leaves,minimum=minimum)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--polynomials',type=Path,default=HERE/'polynomials.json')
    parser.add_argument('--j',type=int,choices=[1,5],required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--nt',type=int,default=94)
    args=parser.parse_args()
    if not __debug__: parser.error('assertions must remain enabled')
    data=json.loads(args.polynomials.read_text())
    coefficients=list(map(Q,data[str(args.j)]['coefficients']))
    started=time.time()
    with mp.Pool(args.workers) as pool:
        outputs=pool.map(produce_one,[(r,coefficients) for r in kr.roots(args.j,args.nt)])
    report=dict(j=args.j,a=str(A),coefficients=list(map(str,coefficients)),nt=args.nt,
                roots=[x[0] for x in outputs],leaves=sum(x[1]['leaves'] for x in outputs),
                evaluations=sum(x[1]['evaluations'] for x in outputs),
                minimum_leaf_lower=min(x[1]['minimum'] for x in outputs),seconds=time.time()-started)
    args.out.write_text(json.dumps(report)+'\n')
    print({k:v for k,v in report.items() if k not in ('roots','coefficients')},flush=True)
