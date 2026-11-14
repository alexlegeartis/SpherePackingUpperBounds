"""Retain the post-archive 44D repair checkpoints with read-back hashes."""
import hashlib,json,tarfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
archive=HERE/'dim44_repair_checkpoints.tar.xz';assert not archive.exists()
keep=set();selected=set()
for n in [22,25,32,*range(35,48)]:
    r=json.loads((ROOT/f'verifications/dim{n}/refinement_profile20261114/certificate.json').read_text())
    source=ROOT/r['research_source'];selected.add(source.parent)
    keep.update([source,source.with_name('candidate.json'),source.with_name('cuts.json'),ROOT/r['angular_file']])
files=[]
for directory in HERE.iterdir():
    if not directory.is_dir() or directory.name in ['sources','__pycache__','dim44_rejected_angular']:continue
    files.extend(p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p not in keep)
manifest={str(p.relative_to(HERE)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(files)}
with tarfile.open(archive,'w:xz',preset=6) as stream:
    for p in sorted(files):stream.add(p,arcname=str(p.relative_to(HERE)),recursive=False)
with tarfile.open(archive,'r:xz') as stream:
    assert len(stream.getmembers())==len(manifest)
    for member in stream:assert hashlib.sha256(stream.extractfile(member).read()).hexdigest()==manifest[member.name]['sha256']
for p in files:
    assert hashlib.sha256(p.read_bytes()).hexdigest()==manifest[str(p.relative_to(HERE))]['sha256'];p.unlink()
for directory in sorted((p for p in HERE.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
    if not any(directory.iterdir()):directory.rmdir()
out={'archive':archive.name,'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
     'archive_bytes':archive.stat().st_size,'uncompressed_bytes':sum(z['bytes'] for z in manifest.values()),
     'files':manifest,'all_members_read_back_and_hashed':True,
     'restore_command':'tar -xJf dim44_repair_checkpoints.tar.xz',
     'scope':'Supplement to ARCHIVE_MANIFEST.json; repaired 44D discovery and the superseded sampled angular source'}
(HERE/'ARCHIVE_MANIFEST_repair_dim44.json').write_text(json.dumps(out,indent=2)+'\n')
index=json.loads((HERE/'DISCOVERY_INDEX.json').read_text())
for entry in index['candidates']:entry['selected_for_frozen_package']=(HERE/entry['path']).parent in selected
seen={entry['path'] for entry in index['candidates']}
for file in HERE.glob('*/candidate.json'):
    if str(file.relative_to(HERE)) in seen:continue
    r=json.loads(file.read_text());index['candidates'].append({'path':str(file.relative_to(HERE)),
       'dimension':r.get('dimension'),'profile_a_exact':r.get('a_exact'),'profile_b_exact':r.get('b_exact'),
       'degree':r.get('degree'),'numerical_density_candidate':r.get('density_candidate'),
       'selected_for_frozen_package':file.parent in selected,'candidate_is_not_a_proof':True})
(HERE/'DISCOVERY_INDEX.json').write_text(json.dumps(index,indent=2)+'\n')
print('REPAIR ARCHIVE VERIFIED',len(files),'files,',out['archive_bytes'],'bytes',flush=True)
