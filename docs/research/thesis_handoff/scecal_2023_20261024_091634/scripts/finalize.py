"""Sequential completion of an interrupted handoff; originals remain read-only."""
import pathlib,subprocess,sys,os,datetime
B=pathlib.Path(__file__).absolute().parents[1]
for name in ['build_bundle.py','finish_indices.py','polish_bundle.py','verify_and_archive.py']:
 print(datetime.datetime.now(datetime.timezone.utc).isoformat(),'START',name,flush=True)
 subprocess.run([sys.executable,'-u',str(B/'scripts'/name)],cwd='/',env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},check=True)
 print('DONE',name,flush=True)
