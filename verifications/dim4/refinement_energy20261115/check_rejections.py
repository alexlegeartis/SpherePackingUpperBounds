"""Meaningful failure controls for certificate coverage, PSD credit and signs."""
from pathlib import Path
from fractions import Fraction as Q
import copy,json,tempfile
import boxes,gram,polynomial

HERE=Path(__file__).resolve().parent


def rejects(name,action):
    try:action()
    except (AssertionError,ValueError):
        print('REJECTED:',name,flush=True)
        return
    raise RuntimeError('invalid certificate accepted: '+name)


def main():
    if not __debug__:raise RuntimeError('assertions must be enabled')
    data=json.loads((HERE/'polynomials.json').read_text())
    coefficients=list(map(Q,data['1']['coefficients']))
    original=json.loads((HERE/'kernel_box_j1.json').read_text())
    with tempfile.TemporaryDirectory(prefix='d4-rejection-controls-') as temp:
        root=Path(temp)
        path=root/'bad.json'
        missing=copy.deepcopy(original);missing['roots'].pop()
        path.write_text(json.dumps(missing))
        rejects('missing root leaves a coverage hole',lambda:boxes.check_tree(path,1,coefficients,1))
        altered=copy.deepcopy(original);altered['roots'][0]['r'][0]='1000001/1000000'
        path.write_text(json.dumps(altered))
        rejects('shifted radial endpoint',lambda:boxes.check_tree(path,1,coefficients,1))
        corrupt=copy.deepcopy(original);corrupt['roots'][0]['tree']=[3,'L','L']
        rejects('invalid subdivision dimension',lambda:boxes.replay_entry((corrupt['roots'][0],coefficients)))
        credit=json.loads((HERE/'gram_credits.json').read_text())
        credit['1'][0]['solution'][0]=str(Q(credit['1'][0]['solution'][0])+1)
        path.write_text(json.dumps(credit))
        rejects('incorrect exact Gram range solution',lambda:gram.check_credits(HERE/'vendor/D4/data/v7/kernel.json',path,1,coefficients))
        rejects('negative coefficients without compensating Gram credits',
                lambda:gram.check_credits(HERE/'vendor/D4/data/v7/kernel.json',None,1,coefficients))
        # A positive constant violates the central sign immediately; this is
        # a whole-domain polynomial check, not a sampled negative control.
        rejects('positive polynomial on the central interval',
                lambda:polynomial.check_nonpositive([Q(1)],Q(133,250),Q(1,2)))
    print('PASS: all six invalid certificates rejected',flush=True)


if __name__=='__main__':main()
