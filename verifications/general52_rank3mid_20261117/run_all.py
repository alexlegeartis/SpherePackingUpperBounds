from pathlib import Path
import hashlib,json,subprocess,sys,tarfile,tempfile
if not __debug__:raise RuntimeError("Assertions must be enabled")
p=Path(__file__).resolve().parent
for name,sha in json.loads((p/'MANIFEST.json').read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==sha,name
with tempfile.TemporaryDirectory(prefix='sphere51-proof-')as tmp:
 root=Path(tmp)
 with tarfile.open(p/'proof.tar.xz','r:xz')as t:
  for member in t.getmembers():
   path=Path(member.name)
   assert not path.is_absolute() and '..'not in path.parts and member.isfile(),member.name
  t.extractall(root,filter='data')
 expected=json.loads((p/'STORED_FILES.json').read_text())
 for name,sha in expected.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name
 subprocess.run([sys.executable,str(root/'verify.py')]+(['--full']if '--full'in sys.argv else []),check=True)
print('PASS portable 52D proof archive'+(' with complete interval replay'if '--full'in sys.argv else ''))
