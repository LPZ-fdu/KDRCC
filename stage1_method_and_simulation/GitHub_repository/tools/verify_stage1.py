"""Verify source/report files and optional GitHub Release assets."""
from pathlib import Path
import argparse,csv,hashlib

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda:f.read(8*1024*1024),b''):h.update(data)
    return h.hexdigest()

def verify(root,assets=None):
    manifests=[(root/'SOURCE_CHECKSUMS.csv',root)]
    if assets is not None:manifests.append((root/'archive/release_assets.csv',Path(assets)))
    count=0
    for manifest,base in manifests:
        with manifest.open(encoding='utf-8-sig',newline='') as f:
            for row in csv.DictReader(f):
                name=row.get('file');p=base/name
                if not p.is_file() or p.stat().st_size!=int(row['bytes']) or digest(p)!=row['sha256']:
                    raise SystemExit('Verification failed: '+str(p))
                count+=1
        print('Verified '+manifest.name,flush=True)
    print('Verified %d files.'%count)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--assets',type=Path);args=a.parse_args()
    verify(Path(__file__).resolve().parents[1],args.assets)
