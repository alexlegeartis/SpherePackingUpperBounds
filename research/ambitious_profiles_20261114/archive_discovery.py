"""Compress retained checkpoints while keeping selected exports easy to replay."""
import hashlib,json,tarfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
assert not (HERE/'discovery_checkpoints.tar.xz').exists(),'Archive already exists'
dimensions=[22,25,32,*range(35,48)]
keep=set();selected=set()
for n in dimensions:
    r=json.loads((ROOT/f'verifications/dim{n}/refinement_profile20261114/certificate.json').read_text())
    source=ROOT/r['research_source'];selected.add(source.parent)
    keep.update([source,source.with_name('candidate.json'),source.with_name('cuts.json'),ROOT/r['angular_file']])
    assert all(p.exists() for p in [source,source.with_name('candidate.json'),ROOT/r['angular_file']])
index=[]
for file in sorted(HERE.glob('*/candidate.json')):
    r=json.loads(file.read_text())
    index.append({'path':str(file.relative_to(HERE)),
                  'dimension':r.get('dimension'),'profile_a_exact':r.get('a_exact'),
                  'profile_b_exact':r.get('b_exact'),'degree':r.get('degree'),
                  'numerical_density_candidate':r.get('density_candidate',r.get('density')),
                  'selected_for_frozen_package':file.parent in selected,
                  'candidate_is_not_a_proof':True})
(HERE/'DISCOVERY_INDEX.json').write_text(json.dumps({'candidates':index,'scope':'Sampled optimization objectives are provisional; consult frozen package replays for certified bounds'},indent=2)+'\n')
files=[]
for directory in HERE.iterdir():
    if not directory.is_dir() or directory.name in ['sources','__pycache__']:continue
    files.extend(p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p not in keep)
manifest={str(p.relative_to(HERE)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(files)}
archive=HERE/'discovery_checkpoints.tar.xz'
with tarfile.open(archive,'w:xz',preset=6) as stream:
    for p in sorted(files):stream.add(p,arcname=str(p.relative_to(HERE)),recursive=False)
with tarfile.open(archive,'r:xz') as stream:
    assert len(stream.getmembers())==len(files)
    for member in stream:
        assert member.isfile() and member.name in manifest
        data=stream.extractfile(member).read();assert hashlib.sha256(data).hexdigest()==manifest[member.name]['sha256']
for p in files:
    assert hashlib.sha256(p.read_bytes()).hexdigest()==manifest[str(p.relative_to(HERE))]['sha256']
    p.unlink()
for directory in sorted((p for p in HERE.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
    if not any(directory.iterdir()):directory.rmdir()
out={'archive':archive.name,'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
     'archive_bytes':archive.stat().st_size,'uncompressed_bytes':sum(z['bytes'] for z in manifest.values()),
     'files':manifest,'all_members_read_back_and_hashed':True,
     'restore_command':'tar -xJf discovery_checkpoints.tar.xz',
     'selected_candidate_coordinates_cuts_and_radial_exports_remain_uncompressed':True}
(HERE/'ARCHIVE_MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')
print('ARCHIVE VERIFIED',len(files),'files,',out['uncompressed_bytes'],'->',out['archive_bytes'],'bytes')
