"""Verify the repository payload against its SHA-256 manifest."""
from pathlib import Path
import argparse
import csv
import hashlib

def verify(manifest,root):
    count=0
    with manifest.open(encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            p=root/row['file']
            if not p.is_file() or p.stat().st_size!=int(row['bytes']):
                raise SystemExit('Missing file or size mismatch: '+str(p))
            h=hashlib.sha256()
            with p.open('rb') as source:
                for block in iter(lambda:source.read(4*1024*1024),b''):
                    h.update(block)
            if h.hexdigest()!=row['sha256']:
                raise SystemExit('SHA-256 mismatch: '+str(p))
            count+=1
    print('Verified',count,'files from',manifest.name)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-assets',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    verify(root/'SUBMISSION_CHECKSUMS.csv',root)
    supplemental_reports=root/'stage1_method_and_simulation/variable_identification_supplement/results'
    if (supplemental_reports/'index.html').is_file():
        verify(supplemental_reports/'results_manifest.csv',supplemental_reports)
    if args.release_assets:
        verify(root/'stage1_method_and_simulation/GitHub_repository/archive/release_assets.csv',args.release_assets)
        supplemental_assets=args.release_assets/'variable_identification_supplement'
        if not (supplemental_assets/'release_assets.csv').is_file():
            supplemental_assets=args.release_assets
        if (supplemental_reports/'index.html').is_file():
            # The supplemental manifest has a unique directory in the local release bundle.
            verify(supplemental_assets/'release_assets.csv',supplemental_assets)
