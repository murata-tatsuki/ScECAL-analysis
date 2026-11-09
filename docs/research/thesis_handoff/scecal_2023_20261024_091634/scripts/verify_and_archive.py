#!/usr/bin/env python3
"""Validate a finished snapshot, hash every bundled file, and archive without symlink traversal."""
import pathlib,json,csv,hashlib,tarfile,os,re,datetime
B=pathlib.Path(__file__).absolute().parents[1]
HASH_CACHE={}
def sha(p):
 key=str(p)
 if key in HASH_CACHE:return HASH_CACHE[key]
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''):h.update(c)
 digest=h.hexdigest();HASH_CACHE[key]=digest
 return digest
def rows(n):
 with (B/n).open() as f:return list(csv.DictReader(f))
m=json.loads((B/'manifest.yaml').read_text())
# Native copies of small calibration ROOTs supplement the previously exported trees.
# Their later capture time remains distinct from the original tree readout time.
native=[]
for record in json.loads((B/'audit/root_readout.json').read_text()):
 if record.get('kind')!='constants' or record.get('size',10**12)>5_000_000:continue
 source_path=pathlib.Path(record['path'])
 try:
  before=source_path.stat();old=datetime.datetime.fromisoformat(record['mtime']).timestamp()
  if before.st_size!=record['size'] or abs(before.st_mtime-old)>0.000002:
   native.append({'path':str(source_path),'availability':'not_copied','reason':'original metadata changed since saved tree export'});continue
  data=source_path.read_bytes();after=source_path.stat()
  if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
   native.append({'path':str(source_path),'availability':'not_copied','reason':'changed during native ROOT read'});continue
  sid='src_'+hashlib.sha256(str(source_path).encode()).hexdigest()[:12]
  relative='evidence/constants/'+sid+'_'+source_path.name;target=B/relative
  target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
  digest=hashlib.sha256(data).hexdigest();captured=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()
  item=next(x for x in m['sources'] if x['source_id']==sid)
  item.update(archive_path=relative,sha256=digest,archive_sha256=digest,copy_sha256=digest,extent='whole_file',availability='copied',captured_at=captured,size=after.st_size,mtime=datetime.datetime.fromtimestamp(after.st_mtime,datetime.timezone(datetime.timedelta(hours=9))).isoformat(),changed_during_read=False,tree_export_captured_at=record.get('read_started_at'),native_copy_note='Later whole-file capture; size/mtime agree with earlier tree export. Snapshot is not atomic.')
  for key in ['hash_reason','hash_omission_reason']:item.pop(key,None)
  native.append({'source_id':sid,'path':str(source_path),'archive_path':relative,'sha256':digest,'captured_at':captured,'availability':'copied'})
 except OSError as exc:native.append({'path':str(source_path),'availability':'not_copied','reason':str(exc)})
m['native_constant_copies']=native
m['methods']['source_epoch']='Most source/tree evidence captured2026-10-24; assembly, verification and selected small native ROOT copies2026-10-26; no atomic or latest-source claim.'
(B/'manifest.yaml').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n')
rr=rows('results.csv');rn=rows('runs.csv');ff=rows('figures.csv')
errors=[];warnings=[]
required=['HANDOFF.md','manifest.yaml','results.csv','runs.csv','figures.csv']
for n in required:
 if not (B/n).is_file():errors.append('missing '+n)
ids={}
for name,key in [('sources','source_id'),('samples','sample_id'),('conditions','condition_id'),('campaigns','campaign_id')]:
 vals=[x[key] for x in m[name]]
 if len(set(vals))!=len(vals):errors.append('duplicate '+name+' IDs')
 ids[name]=set(vals)
for table in [rr,rn,ff]:
 for r in table:
  for key,collection in [('campaign_id','campaigns'),('sample_id','samples'),('condition_id','conditions')]:
   v=r.get(key,'')
   if v and v!='unknown' and v not in ids[collection]:errors.append('unresolved '+key+' '+v)
for r in rr:
 if r['valid'] not in ['true','false','unknown']:errors.append('bad valid '+r['result_id'])
 if r['status']!='暫定':warnings.append('non-default adoption '+r['result_id'])
 if r['valid']=='false' and not r['invalid_reason']:errors.append('missing invalid reason '+r['result_id'])
 if r['value']=='' and not r['invalid_reason']:errors.append('empty value without explanation '+r['result_id'])
 if re.search(r'/202[25]/',r['source']):errors.append('out-of-scope measurement result '+r['source'])
for c in m['conditions']:
 for sid in c.get('sample_ids',[]):
  if sid!='unknown' and sid not in ids['samples']:errors.append('condition sample missing '+sid)
 for sid in c.get('source_ids',[]):
  if sid!='unknown' and sid not in ids['sources']:errors.append('condition source missing '+sid)
for s in m['sources']:
 rel=s.get('archive_path')
 if not rel:continue
 p=B/rel
 if not p.is_file():errors.append('missing source copy '+rel);continue
 h=sha(p);exp=s.get('archive_sha256',s.get('sha256'))
 if exp and h!=exp:errors.append('source copy hash mismatch '+rel)
 # Compare original/copy hash recorded at capture without rereading mutable originals.
 old=s.get('sha256');cop=s.get('archive_sha256')
 if old and cop and old!=cop and s.get('extent')=='whole_file':errors.append('capture original/copy mismatch '+rel)
for p in B.rglob('*'):
 if p.is_symlink():errors.append('symlink in snapshot '+str(p.relative_to(B)))
 if p.is_file() and p.stat().st_size>100_000_000:errors.append('unexpected large file '+str(p.relative_to(B)))
for f in ff:
 if f.get('availability')=='copied':
  p=B/f['archive_path']
  if not p.is_file():errors.append('missing figure '+str(p))
  elif f.get('sha256') and sha(p)!=f['sha256']:errors.append('figure hash mismatch '+f['figure_id'])
text=(B/'HANDOFF.md').read_text()+(B/'SOURCES.md').read_text()
for target in re.findall(r'\]\(([^)]+)\)',text):
 if target.startswith(('http:','https:','#')):continue
 if not (B/target.split('#')[0]).exists():errors.append('broken relative link '+target)
report={'validated_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest_parser':'JSON-compatible YAML 1.2','results':len(rr),'runs':len(rn),'figures':len(ff),'errors':errors,'warnings':warnings,'new_fits_or_event_processing':False,'original_files_modified':False,'scope_check':'2022/2025 measurements excluded; cross-year reference documents separately retained','validation':'internal IDs, missing/invalid handling, copied-file hashes, relative links, no symlinks/large raw files'}
(B/'audit/validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
if errors:print(json.dumps(report,ensure_ascii=False,indent=2));raise SystemExit(1)
# All writers must have completed before this step. The manifest itself is covered by this inventory.
files=sorted(p for p in B.rglob('*') if p.is_file() and p != B/'SHA256SUMS.txt')
hashes={str(p.relative_to(B)):sha(p) for p in files}
(B/'SHA256SUMS.txt').write_text(''.join(h+'  '+rel+'\n' for rel,h in hashes.items()))
archive=B.parent/(B.name+'.tar.gz')
if archive.exists():raise RuntimeError('refusing to overwrite existing final archive')
partial=archive.parent/(archive.name+'.part')
if partial.exists():partial.unlink()
with tarfile.open(partial,'w:gz',compresslevel=6,dereference=False) as t:t.add(B,arcname=B.name,recursive=True)
count=0
with tarfile.open(partial,'r:gz') as t:
 for x in t:
  if x.issym() or x.islnk():raise RuntimeError('archive contains link '+x.name)
  if not x.isfile():continue
  rel=str(pathlib.PurePosixPath(x.name).relative_to(B.name))
  if rel=='SHA256SUMS.txt':continue
  h=hashlib.sha256(t.extractfile(x).read()).hexdigest()
  if hashes.get(rel)!=h:raise RuntimeError('archive content hash mismatch '+rel)
  count+=1
if count!=len(hashes):raise RuntimeError('archive inventory mismatch')
archive_sha=sha(partial)
partial.replace(archive)
(archive.parent/(archive.name+'.sha256')).write_text(archive_sha+'  '+archive.name+'\n')
print(json.dumps({'snapshot_path':str(B),'snapshot_file_bytes':sum(p.stat().st_size for p in B.rglob('*') if p.is_file()),'archive_path':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':archive_sha,'archive_checksum_path':str(archive)+'.sha256','hashed_files':len(hashes),'tar_contents_verified':True,'results':len(rr),'runs':len(rn),'figures':len(ff)},ensure_ascii=False,indent=2),flush=True)
