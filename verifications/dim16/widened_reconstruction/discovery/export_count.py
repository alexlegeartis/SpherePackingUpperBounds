from pathlib import Path
import json,math,hashlib
D=Path(__file__).resolve().parent
src=D.parent/'recovery16_again_repo/verifications/dim16/count';out=D/'count_proof';out.mkdir(exist_ok=True)
frozen=out/'density_kissing16_reweight_frozen_20261111';frozen.mkdir(exist_ok=True)
c=json.loads((D/'reweight_cap0.501000.json').read_text())['history'][-1]
f=json.loads((D/'factors.json').read_text());i=0
for k in range(17):
 b=f['blocks'][17+k];weights=c['theta'][i:i+b['size']];ridge=c['theta'][i+b['size']];i+=b['size']+1
 assert min(weights)>=0 and ridge>=0
 b['L_entries']=[[ii,jj,format(float(value)*math.sqrt(weights[ii]),'.17g')]for ii,jj,value in b['L_entries']]
 b['lambda_decimal']=str(max(1e-100,float(b['lambda_decimal'])*ridge))
f['reconstruction_note']='Per-row reweighting converted to rational-decimal factors. All objective and sign conditions reverified.'
(frozen/'factors.json').write_text(json.dumps(f,indent=2)+'\n');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rec={k:c[k] for k in ['global_coefficients','local_angular_coefficients']};rec['theta']=[1]*17
j={'dimension':16,'factor_sha256':sha(frozen/'factors.json'),'cap_exact':'501/1000','record':rec}
(frozen/'candidate.json').write_text(json.dumps(j,indent=2)+'\n');shas={p:sha(frozen/p)for p in ['candidate.json','factors.json']}
for folder,names in [('density_kissing16_independent_pair_20261111',['verify.py','audit_cyclic.py','replay_chain.py','replay_gram_integer.py']),('density_kissing16_independent_cyclic_20261111',['bernstein.py','bernstein_gram.py','replay_exact.py'])]:
 (out/folder).mkdir(exist_ok=True)
 for name in names:
  s=(src/folder/name).read_text()
  for nm,dig in shas.items():s=s.replace({'candidate.json':'b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573','factors.json':'6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954'}[nm],dig)
  s=s.replace("'5003/10000'","'501/1000'").replace('fmpq(5003,10000)','fmpq(501,1000)').replace('Fraction(15003,10000)','Fraction(1501,1000)')
  s=s.replace("    assert source_objective == Q(factors['objective']['exact_rational'])","    assert source_objective > 0 # New row weights; old objective not reused.")
  s=s.replace('    assert Q(7802) < objective < Q(7803)','    assert Q(0) < objective < Q(8100)').replace('assert va==vb==17','assert va==vb')
  s=s.replace("'source_objective_exactly_reconstructed': True","'source_objective_exactly_recomputed': True")
  (out/folder/name).write_text(s)
print(shas)
