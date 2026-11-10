import sys,json
from fractions import Fraction as F
sys.path.insert(0,'solver_packages')
from flint import arb,fmpq,ctx
ctx.prec=256;p=sys.argv[1];a=json.load(open(p));b=json.load(open(p.replace('.json','_psd_exact.json')))
c=[fmpq(str(F(str(z))))for z in a['x'][:32]];h=sum(((-1)**i*v for i,v in enumerate(c)),fmpq(0));num=sum(c,fmpq(0))+fmpq(b['matrices'][0][0][0])+sum((fmpq(str(F(str(w))))*fmpq(str(F(str(v))))for w,v in zip(a['count_weights'],a['count_bounds'])),fmpq(0));v=arb.pi()**2/60*arb(num)/arb(h);out={'density_enclosure':str(v),'below_0.510584101':bool(v<arb('0.510584101')),'normalization_exact':str(h),'numerator_exact':str(num),'all_count_weights_nonnegative':all(w>=0 for w in a['count_weights'])};print(json.dumps(out));json.dump(out,open(p.replace('.json','_density_exact.json'),'w'),indent=2)
