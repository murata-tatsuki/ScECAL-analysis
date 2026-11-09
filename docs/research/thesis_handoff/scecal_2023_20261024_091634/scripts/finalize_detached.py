"""Resume this handoff across chat restarts; never alter original analysis files."""
import pathlib,subprocess,sys,os,json,time,datetime
B=pathlib.Path(__file__).absolute().parents[1]
job={'pid':os.getpid(),'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wait_for_prior_worker':{'pid':1059848,'start_ticks':'88293161'},'scope':'Only this unfinished handoff, its archive and checksum','stages':['build_bundle.py','finish_indices.py','polish_bundle.py','verify_and_archive.py']}
(B/'audit/finalize_job.json').write_text(json.dumps(job,indent=2)+'\n')
while True:
 try:
  state=pathlib.Path('/proc/1059848/stat').read_text().split()
  if state[21]!='88293161' or state[2]=='Z':break
 except FileNotFoundError:break
 time.sleep(5)
try:
 with (B/'audit/finalize_run.log').open('w') as log:
  for name in job['stages'][:3]:
   print(datetime.datetime.now(datetime.timezone.utc).isoformat(),'START',name,file=log,flush=True)
   subprocess.run([sys.executable,'-u',str(B/'scripts'/name)],cwd='/',env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},stdout=log,stderr=subprocess.STDOUT,check=True)
   print('DONE',name,file=log,flush=True)
  print('Preparation complete; log closed before checksumming and archiving.',file=log,flush=True)
 # No snapshot file may be appended after the verifier creates SHA256SUMS.
 r=subprocess.run([sys.executable,'-u',str(B/'scripts/verify_and_archive.py')],cwd='/',env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 if r.returncode:
  (B/'audit/finalize_failure.log').write_text(r.stdout);raise SystemExit(r.returncode)
except Exception as e:
 (B/'audit/finalize_failure.log').write_text(type(e).__name__+': '+str(e)+'\n');raise
