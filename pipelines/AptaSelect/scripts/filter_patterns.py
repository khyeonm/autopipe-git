#!/usr/bin/env python3
"""Stages 2-4: Nested pattern filtering with fuzzy matching."""

import sys
import argparse

def hamming(s1, s2):
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_pattern(seq, pattern, max_mismatches, start=0):
    plen = len(pattern)
    for i in range(start, len(seq) - plen + 1):
        if hamming(seq[i:i+plen], pattern) <= max_mismatches:
            return i
    return -1

def filter_sequences(input_path, output_path, pat1, pat2, max_mm, between_length=None):
    counts = {}
    total = 0
    passed = 0

    with open(input_path) as f:
        for line in f:
            seq = line.strip()
            if not seq:
                continue
            total += 1

            pos1 = find_pattern(seq, pat1, max_mm)
            if pos1 < 0:
                continue

            pos2 = find_pattern(seq, pat2, max_mm, start=pos1 + len(pat1))
            if pos2 < 0:
                continue

            # Check between-length if required
            if between_length is not None:
                between = pos2 - (pos1 + len(pat1))
                if between != between_length:
                    continue

            # Pattern-inclusive extraction
            extracted = seq[pos1: pos2 + len(pat2)]
            counts[extracted] = counts.get(extracted, 0) + 1
            passed += 1

    print(f"Total: {total}, Passed: {passed}, Unique: {len(counts)}", file=sys.stderr)

    with open(output_path, 'w') as out:
        for seq, count in counts.items():
            out.write(f"{seq}\t{count}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--pat1', required=True)
    parser.add_argument('--pat2', required=True)
    parser.add_argument('--max-mismatches', type=int, default=1)
    parser.add_argument('--between-length', type=int, default=None)
    parser.add_argument('--input-has-counts', action='store_true',
                        help='Input is tab-separated seq\\tcount (from previous stage)')
    args = parser.parse_args()

    counts = {}
    total = 0
    passed = 0

    with open(args.input) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if args.input_has_counts:
                parts = line.split('\t')
                seq = parts[0]
                weight = int(parts[1]) if len(parts) > 1 else 1
            else:
                seq = line
                weight = 1

            total += weight

            pos1 = find_pattern(seq, args.pat1, args.max_mismatches)
            if pos1 < 0:
                continue

            pos2 = find_pattern(seq, args.pat2, args.max_mismatches, start=pos1 + len(args.pat1))
            if pos2 < 0:
                continue

            if args.between_length is not None:
                between = pos2 - (pos1 + len(args.pat1))
                if between != args.between_length:
                    continue

            extracted = seq[pos1: pos2 + len(args.pat2)]
            counts[extracted] = counts.get(extracted, 0) + weight
            passed += weight

    print(f"Total: {total}, Passed: {passed}, Unique: {len(counts)}", file=sys.stderr)

    with open(args.output, 'w') as out:
        for seq, count in counts.items():
            out.write(f"{seq}\t{count}\n")

if __name__ == '__main__':
    main()