#!/usr/bin/env python3
"""Read-only source capture. Writes only beneath this new snapshot."""
import argparse, concurrent.futures, hashlib, json, os, pathlib, shutil
from datetime import datetime
from zoneinfo import ZoneInfo
BASE = pathlib.Path(__file__).resolve().parents[1]
def now(): return datetime.now(ZoneInfo('Asia/Tokyo')).isoformat()
def capture(item):
    p = pathlib.Path(item['path'])
    rec = dict(item, captured_at=now(), evidence_status='確認済み')
    rec['source_id'] = 'src_' + hashlib.sha256(str(p).encode()).hexdigest()[:12]
    try:
        a = p.stat(); rec.update(size=a.st_size,mtime=datetime.fromtimestamp(a.st_mtime,ZoneInfo('Asia/Tokyo')).isoformat(),resolved_path=os.path.realpath(p))
        if a.st_size > item.get('max_size',5000000):
            rec.update(availability='metadata_only',sha256=None,not_copied_reason='size_limit; no full raw/event ROOT copy')
            return rec
        data = p.read_bytes(); rec['sha256']=hashlib.sha256(data).hexdigest()
        rel = item.get('copy_to') or ('evidence/captured/'+rec['source_id']+'_'+p.name)
        target=BASE/rel; target.parent.mkdir(parents=True,exist_ok=True)
        if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest()!=rec['sha256']:target.write_bytes(data)
        z=p.stat(); rec.update(archive_path=rel,availability='copied',capture_extent='whole',changed_during_read=(a.st_size!=z.st_size or a.st_mtime_ns!=z.st_mtime_ns),finished_at=now(),copy_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    except Exception as e: rec.update(availability='unavailable',error=str(e),sha256=None)
    return rec
if __name__=='__main__':
    q=argparse.ArgumentParser(); q.add_argument('plan'); q.add_argument('output'); args=q.parse_args()
    items=json.loads((BASE/args.plan).read_text()); out=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for r in pool.map(capture,items):
            out.append(r)
            if len(out)%25==0: print('captured',len(out),flush=True)
    (BASE/args.output).write_text(json.dumps(out,ensure_ascii=False,indent=2))
    print('done',len(out),sum(r.get('size',0) for r in out if r['availability']=='copied'),flush=True)
