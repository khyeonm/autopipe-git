#!/usr/bin/env python3
"""Split a joined-reads file into chunks for parallel processing."""

import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--chunk-size', type=int, default=10000)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    chunk_idx = 0
    lines = []

    with open(args.input) as f:
        for line in f:
            lines.append(line)
            if len(lines) >= args.chunk_size:
                out_path = os.path.join(args.output_dir, f"chunk_{chunk_idx:06d}.txt")
                with open(out_path, 'w') as out:
                    out.writelines(lines)
                lines = []
                chunk_idx += 1

    if lines:
        out_path = os.path.join(args.output_dir, f"chunk_{chunk_idx:06d}.txt")
        with open(out_path, 'w') as out:
            out.writelines(lines)
        chunk_idx += 1

    print(f"Created {chunk_idx} chunks", file=sys.stderr)

    # Write chunk list for Snakemake
    manifest = os.path.join(args.output_dir, "chunks.txt")
    with open(manifest, 'w') as m:
        for i in range(chunk_idx):
            m.write(f"chunk_{i:06d}\n")

if __name__ == '__main__':
    main()