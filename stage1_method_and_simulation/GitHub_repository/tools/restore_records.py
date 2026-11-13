"""Verify and restore independently readable stage-one ZIP assets."""
from pathlib import Path,PurePosixPath
import argparse,csv,hashlib,zipfile,shutil
from verify_stage1 import digest

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--assets',required=True,type=Path)
    a.add_argument('--destination',required=True,type=Path);args=a.parse_args()
    root=Path(__file__).resolve().parents[1];dest=args.destination.resolve()
    dest.mkdir(parents=True,exist_ok=True)
    with (root/'archive/release_assets.csv').open(encoding='utf-8-sig',newline='') as f:assets=list(csv.DictReader(f))
    for row in assets:
        path=args.assets/row['file']
        if path.stat().st_size!=int(row['bytes']) or digest(path)!=row['sha256']:
            raise SystemExit('Asset verification failed: '+str(path))
        with zipfile.ZipFile(path) as z:
            for info in z.infolist():
                parts=PurePosixPath(info.filename).parts
                if not parts or parts[0]!='KDRCC_formal_20261027' or '..' in parts or '\\' in info.filename:
                    raise SystemExit('Invalid archive path: '+info.filename)
                target=dest.joinpath(*parts).resolve()
                if dest not in target.parents:raise SystemExit('Path outside destination: '+str(target))
                if target.exists():
                    with z.open(info) as src:
                        h=hashlib.sha256()
                        for data in iter(lambda:src.read(1024*1024),b''):h.update(data)
                    if digest(target)!=h.hexdigest():raise SystemExit('Existing file differs: '+str(target))
                else:
                    target.parent.mkdir(parents=True,exist_ok=True)
                    with z.open(info) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst,1024*1024)
        print('Restored '+row['file'],flush=True)
    print('Records restored to '+str(dest/'KDRCC_formal_20261027'))
