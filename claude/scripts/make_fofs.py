#!/usr/bin/env python3
"""Split a list of sample files into C chunk FOFs in kmtricks format ("id : path").

  ./make_fofs.py samples.txt -c 4 -o fofs        # one path per line
  ./make_fofs.py /data/unitigs -c 4 -o fofs      # directory scanned recursively
  ./make_fofs.py samples.txt -c 1 -n 250 -o fofs # keep only the first 250 samples

Sample ids are file basenames without extensions (unique across chunks is required
by kmindex merge). Chunks are contiguous and balanced in count, not in size.
"""

import argparse
import os
import re

EXT = re.compile(r"\.(fa|fasta|fna|fq|fastq)(\.(gz|zst|bz2))?$")


def collect(src):
    if os.path.isdir(src):
        paths = []
        for root, _, files in os.walk(src):
            paths += [os.path.join(root, f) for f in files if EXT.search(f)]
        return sorted(paths)
    with open(src) as f:
        return [l.strip() for l in f if l.strip() and not l.startswith("#")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", help="directory or text file listing sample files")
    ap.add_argument("-c", "--chunks", type=int, default=1)
    ap.add_argument("-n", "--max-samples", type=int, default=0)
    ap.add_argument("-o", "--output", required=True, help="output directory")
    a = ap.parse_args()

    paths = collect(a.source)
    if a.max_samples:
        paths = paths[: a.max_samples]
    if not paths:
        raise SystemExit("no sample found")
    ids = [EXT.sub("", os.path.basename(p)) for p in paths]
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate sample ids")

    os.makedirs(a.output, exist_ok=True)
    size = -(-len(paths) // a.chunks)
    for c in range(a.chunks):
        block = list(zip(ids, paths))[c * size:(c + 1) * size]
        if not block:
            break
        with open(os.path.join(a.output, f"chunk_{c}.fof"), "w") as f:
            for sid, p in block:
                f.write(f"{sid} : {os.path.abspath(p)}\n")
        print(f"chunk_{c}.fof: {len(block)} samples")


if __name__ == "__main__":
    main()
