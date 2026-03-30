#!/usr/bin/env python3
"""
AptaSelect Stage 1: Join paired-end reads.
Reverse-complements Read 2, finds best overlap, merges reads using quality scores.
"""

import sys
import json
from itertools import zip_longest

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def reverse_complement(seq):
    return seq.translate(COMP)[::-1]

def parse_fastq_pairs(r1_path, r2_path):
    """Yield (name, seq1, qual1, seq2, qual2) tuples from paired FASTQ files."""
    with open(r1_path) as f1, open(r2_path) as f2:
        while True:
            lines1 = [f1.readline().rstrip('\n') for _ in range(4)]
            lines2 = [f2.readline().rstrip('\n') for _ in range(4)]
            if not lines1[0] or not lines2[0]:
                break
            yield (lines1[0], lines1[1], lines1[3], lines2[1], lines2[3])

def find_best_overlap(rd1, qual1, rd2, qual2, min_overlap, max_pct_diff):
    """Find the best overlap between rd1's tail and rd2's head."""
    max_ov = min(len(rd1), len(rd2))
    best_score = float('inf')
    best_ov = -1
    best_d = 0

    for ov in range(min_overlap, max_ov + 1):
        # Compare tail of rd1 with head of rd2
        seg1 = rd1[len(rd1) - ov:]
        seg2 = rd2[:ov]
        d = sum(1 for a, b in zip(seg1, seg2) if a != b)
        if d / ov <= max_pct_diff / 100.0:
            score = (1000 * (d * d + 1)) / ov
            if score < best_score:
                best_score = score
                best_ov = ov
                best_d = d

    return best_ov, best_d

def merge_reads(rd1, qual1, rd2, qual2, ov):
    """Merge two reads using overlap, resolving mismatches by quality."""
    if ov <= 0:
        return rd1 + rd2

    overlap_start = len(rd1) - ov
    merged = list(rd1)

    # Resolve overlap mismatches using quality scores
    for i in range(ov):
        pos1 = overlap_start + i
        pos2 = i
        if rd1[pos1] != rd2[pos2]:
            q1 = ord(qual1[pos1]) - 33
            q2 = ord(qual2[pos2]) - 33
            if q2 > q1:
                merged[pos1] = rd2[pos2]

    # Append non-overlapping part of rd2
    merged_seq = ''.join(merged) + rd2[ov:]
    return merged_seq

def main():
    r1_path = sys.argv[1]
    r2_path = sys.argv[2]
    output_path = sys.argv[3]
    config_path = sys.argv[4]

    with open(config_path) as f:
        config = json.load(f)

    min_overlap = config['min_overlap']
    pct_diff = config['pct_diff']
    is_short = config['is_short']

    total = 0
    joined = 0

    with open(output_path, 'w') as out:
        for name, seq1, qual1, seq2, qual2 in parse_fastq_pairs(r1_path, r2_path):
            total += 1

            rc_seq2 = reverse_complement(seq2)
            rc_qual2 = qual2[::-1]

            if is_short:
                # Short library: rd1 = RC(Read2), rd2 = Read1
                rd1, q1 = rc_seq2, rc_qual2
                rd2, q2 = seq1, qual1
            else:
                # Long library: rd1 = Read1, rd2 = RC(Read2)
                rd1, q1 = seq1, qual1
                rd2, q2 = rc_seq2, rc_qual2

            ov, d = find_best_overlap(rd1, q1, rd2, q2, min_overlap, pct_diff / 100.0 * 100)
            # Re-check: the find_best_overlap already uses max_pct_diff correctly
            if ov < 0:
                continue

            merged = merge_reads(rd1, q1, rd2, q2, ov)
            out.write(merged + '\n')
            joined += 1

    # Write stats
    stats_path = output_path + '.stats'
    with open(stats_path, 'w') as sf:
        sf.write(f"total_pairs\t{total}\n")
        sf.write(f"joined\t{joined}\n")
        sf.write(f"join_rate\t{joined/total*100:.2f}%\n" if total > 0 else "join_rate\t0%\n")

    print(f"Join: {joined}/{total} pairs joined ({joined/total*100:.2f}%)" if total > 0 else "Join: 0 pairs")

if __name__ == '__main__':
    main()