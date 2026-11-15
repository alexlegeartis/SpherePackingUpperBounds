"""Exercise modified historical registry/hash gates without repeating domains."""
from pathlib import Path
import ast,json,hashlib

if not __debug__:raise RuntimeError('assertions must be enabled')
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
names=['verifications/dim17/verify.py','verifications/dim17/refinement_repaired/verify.py',
       'verifications/dim18/verify.py','verifications/dim18/refinement_cyclic/verify.py',
       'verifications/dim19/refinement_cyclic/verify.py','verifications/dim20/refinement_cyclic/verify.py',
       'verifications/dim21/refinement_gaussian/verify.py']
results=[]
for name in names:
    file=ROOT/name;text=file.read_text();tree=ast.parse(text)
    def replay(node):
        return any(isinstance(x,ast.Call) and isinstance(x.func,ast.Attribute) and x.func.attr=='run'
                   and isinstance(x.func.value,ast.Name) and x.func.value.id=='subprocess' for x in ast.walk(node))
    cut=next(node.lineno for node in tree.body if replay(node))
    prefix='\n'.join(text.splitlines()[:cut-1])
    exec(compile(prefix,str(file),'exec'),{'__file__':str(file),'__name__':'__main__'})
    results.append(dict(path=name,pass_prefix=True,sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    print('HISTORICAL PREFIX PASS',name,flush=True)
(HERE/'HISTORICAL_PREFIX_CHECKS.json').write_text(json.dumps(dict(pass_all=True,scope='Registry and hash prefixes; unchanged mathematical domain proofs are not rerun here.',results=results),indent=2)+'\n')
