"""Download Hu et al. single-cell analysis inputs and results with SHA-256 checks."""
import argparse
import csv
import hashlib
from pathlib import Path
import urllib.request

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include', action='append', choices=['data', 'results', 'supplementary_tables'])
    args = parser.parse_args()
    manifest = Path(__file__).resolve().parents[1] / 'config' / 'PUBLIC_DATA_FILES.tsv'
    with manifest.open(newline='') as source:
        rows = list(csv.DictReader(source, delimiter='\t'))
    root = args.output.resolve()
    for row in rows:
        relative = Path(row['path'])
        if args.include and relative.parts[0] not in args.include:
            continue
        target = (root / relative).resolve()
        if root not in target.parents:
            raise ValueError('Invalid manifest path: ' + row['path'])
        if target.is_file() and target.stat().st_size == int(row['bytes']) and sha256(target) == row['sha256']:
            print('Verified existing ' + row['path'], flush=True)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + '.part')
        digest = hashlib.sha256()
        print('Downloading ' + row['path'], flush=True)
        with urllib.request.urlopen(row['url'], timeout=120) as response, partial.open('wb') as destination:
            while True:
                chunk = response.read(8 * 1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
                destination.write(chunk)
        if partial.stat().st_size != int(row['bytes']) or digest.hexdigest() != row['sha256']:
            raise RuntimeError('Download checksum mismatch: ' + row['path'])
        partial.replace(target)
    print('Selected files verified.', flush=True)

if __name__ == '__main__':
    main()
