"""Verify the distributed SHA-256 file inventory (Python standard library)."""
import csv,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
fail=[];count=0
with (root/'SOURCE_CHECKSUMS.csv').open(encoding='utf-8-sig',newline='') as f:
    for row in csv.DictReader(f):
        path=root/row['file'];count+=1
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
            fail.append(row['file'])
if fail:raise SystemExit('Checksum mismatch: '+', '.join(fail))
print('Verified %d release files.'%count)
