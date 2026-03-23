#!/usr/bin/env python3
"""Stage 1: Join paired-end reads using overlap detection."""

import sys
import gzip
import argparse
from itertools import islice

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def phred_to_qual(c):
    return ord(c) - 33

def join_pair(r1_seq, r1_qual, r2_seq, r2_qual, min_overlap, pct_diff, is_short):
    rc_r2 = reverse_complement(r2_seq)
    rc_r2_qual = r2_qual[::-1]

    if is_short:
        rd1_seq, rd1_qual = rc_r2, rc_r2_qual
        rd2_seq, rd2_qual = r1_seq, r1_qual
    else:
        rd1_seq, rd1_qual = r1_seq, r1_qual
        rd2_seq, rd2_qual = rc_r2, rc_r2_qual

    best_score = None
    best_ov = None

    max_ov = min(len(rd1_seq), len(rd2_seq))
    for ov in range(min_overlap, max_ov + 1):
        tail = rd1_seq[-ov:]
        head = rd2_seq[:ov]
        d = hamming(tail, head)
        if d / ov <= pct_diff / 100:
            score = (1000 * (d * d + 1)) / ov
            if best_score is None or score < best_score:
                best_score = score
                best_ov = ov

    if best_ov is None:
        return None

    ov = best_ov
    # Resolve mismatches using quality scores
    merged_overlap = []
    for i in range(ov):
        b1 = rd1_seq[-(ov - i)]
        b2 = rd2_seq[i]
        q1 = phred_to_qual(rd1_qual[-(ov - i)])
        q2 = phred_to_qual(rd2_qual[i])
        merged_overlap.append(b1 if q1 >= q2 else b2)

    joined = rd1_seq[:-ov] + ''.join(merged_overlap) + rd2_seq[ov:]
    return joined

def open_file(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')

def read_fastq(fh):
    while True:
        header = fh.readline().strip()
        if not header:
            break
        seq = fh.readline().strip()
        fh.readline()  # +
        qual = fh.readline().strip()
        yield seq, qual

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-overlap', type=int, default=6)
    parser.add_argument('--pct-diff', type=float, default=8)
    parser.add_argument('--is-short', type=lambda x: x.lower() == 'true', default=True)
    args = parser.parse_args()

    joined = 0
    failed = 0
    with open_file(args.r1) as f1, open_file(args.r2) as f2, open(args.output, 'w') as out:
        for (s1, q1), (s2, q2) in zip(read_fastq(f1), read_fastq(f2)):
            result = join_pair(s1, q1, s2, q2, args.min_overlap, args.pct_diff, args.is_short)
            if result:
                out.write(result + '\n')
                joined += 1
            else:
                failed += 1

    print(f"Joined: {joined}, Failed: {failed}", file=sys.stderr)

if __name__ == '__main__':
    main()