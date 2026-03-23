#!/usr/bin/env python3
"""Stage 5: Aggregate counts across chunks and output top-N sequences."""

import sys
import argparse
import glob
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', required=True, help='Directory with stage4 chunk outputs')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--top-n', type=int, default=10)
    args = parser.parse_args()

    aggregated = {}

    chunk_files = sorted(glob.glob(os.path.join(args.input_dir, 'chunk_*.tsv')))
    if not chunk_files:
        # Try reading a single file
        chunk_files = [args.input_dir] if os.path.isfile(args.input_dir) else []

    for fpath in chunk_files:
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                seq = parts[0]
                count = int(parts[1])
                aggregated[seq] = aggregated.get(seq, 0) + count

    ranked = sorted(aggregated.items(), key=lambda x: x[1], reverse=True)
    top = ranked[:args.top_n]

    with open(args.output, 'w') as out:
        out.write("rank\tsequence\tcount\n")
        for rank, (seq, count) in enumerate(top, 1):
            out.write(f"{rank}\t{seq}\t{count}\n")

    print(f"Total unique sequences: {len(aggregated)}", file=sys.stderr)
    print(f"Top {args.top_n} written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()