#!/usr/bin/env python3
"""
AptaSelect: Identify high-frequency aptamer candidate sequences from paired-end SELEX FASTQ files.
Stages: Join -> Selection Filtering -> 1st Sort Filtering -> 2nd Sort Filtering -> Aggregation & Ranking
"""

import sys
import gzip
import json
import os
from collections import Counter

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def parse_fastq(filepath):
    """Yield (name, seq, qual) tuples from a FASTQ file (gzipped or plain)."""
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            yield header, seq, qual

def phred_score(char):
    return ord(char) - 33

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def join_reads(seq1, qual1, seq2, qual2, min_overlap, pct_diff, is_short):
    """
    Join paired-end reads by finding the best overlap.
    Returns joined sequence or None if no valid overlap found.
    """
    rc_seq2 = reverse_complement(seq2)
    rc_qual2 = qual2[::-1]

    if is_short:
        # Short library: rd1 = RC(Read2), rd2 = Read1
        rd1_seq, rd1_qual = rc_seq2, rc_qual2
        rd2_seq, rd2_qual = seq1, qual1
    else:
        # Long library: rd1 = Read1, rd2 = RC(Read2)
        rd1_seq, rd1_qual = seq1, qual1
        rd2_seq, rd2_qual = rc_seq2, rc_qual2

    max_ov = min(len(rd1_seq), len(rd2_seq))
    best_score = float('inf')
    best_ov = -1

    for ov in range(min_overlap, max_ov + 1):
        tail = rd1_seq[-ov:]
        head = rd2_seq[:ov]
        d = hamming_distance(tail, head)
        if d / ov <= pct_diff / 100.0:
            score = (1000.0 * (d * d + 1)) / ov
            if score < best_score:
                best_score = score
                best_ov = ov

    if best_ov < 0:
        return None

    # Resolve mismatches using quality scores
    ov = best_ov
    merged_overlap = []
    for i in range(ov):
        pos_rd1 = len(rd1_seq) - ov + i
        pos_rd2 = i
        if rd1_seq[pos_rd1] == rd2_seq[pos_rd2]:
            merged_overlap.append(rd1_seq[pos_rd1])
        else:
            if phred_score(rd1_qual[pos_rd1]) >= phred_score(rd2_qual[pos_rd2]):
                merged_overlap.append(rd1_seq[pos_rd1])
            else:
                merged_overlap.append(rd2_seq[pos_rd2])

    joined = rd1_seq[:-ov] + ''.join(merged_overlap) + rd2_seq[ov:]
    return joined

def pattern_matches(seq, pattern, max_mismatches):
    """Yield all start positions where pattern matches seq with <= max_mismatches."""
    plen = len(pattern)
    for i in range(len(seq) - plen + 1):
        if hamming_distance(seq[i:i+plen], pattern) <= max_mismatches:
            yield i

def find_pattern_pair(seq, pattern1, pattern2, max_mismatches, required_between=None):
    """
    Find the first valid (p1_pos, p2_pos) pair using nested loop.
    Returns (p1_start, p2_end) or None.
    p2 must be downstream of p1.
    If required_between is set, the region between patterns must have exactly that length.
    """
    p1len = len(pattern1)
    p2len = len(pattern2)

    for p1_pos in pattern_matches(seq, pattern1, max_mismatches):
        after_p1 = p1_pos + p1len
        for p2_pos in pattern_matches(seq[after_p1:], pattern2, max_mismatches):
            actual_p2_pos = after_p1 + p2_pos
            if required_between is not None:
                between_len = actual_p2_pos - after_p1
                if between_len != required_between:
                    continue
            return p1_pos, actual_p2_pos + p2len
        # If we iterated all p2 positions for this p1 without success, continue to next p1
    return None

def process_chunk(records, config):
    """Process a chunk of paired-end records through all stages."""
    min_overlap = config['min_overlap']
    pct_diff = config['pct_diff']
    is_short = config['is_short']
    sel_read1 = config['sel_read1']
    sel_read2 = config['sel_read2']
    s1_read1 = config['s1_read1']
    s1_read2 = config['s1_read2']
    s1_length = config['s1_length']
    s2_read1 = config['s2_read1']
    s2_read2 = config['s2_read2']
    s2_length = config['s2_length']
    max_mismatches = config['max_mismatches']

    counts = {
        'total': 0,
        'joined': 0,
        'selection': 0,
        'sort1': 0,
        'sort2': 0,
    }
    seq_counter = Counter()

    for seq1, qual1, seq2, qual2 in records:
        counts['total'] += 1

        # Stage 1: Join
        joined = join_reads(seq1, qual1, seq2, qual2, min_overlap, pct_diff, is_short)
        if joined is None:
            continue
        counts['joined'] += 1

        # Stage 2: Selection Filtering
        result = find_pattern_pair(joined, sel_read1, sel_read2, max_mismatches, required_between=None)
        if result is None:
            continue
        p1_start, p2_end = result
        sel_seq = joined[p1_start:p2_end]
        counts['selection'] += 1

        # Stage 3: 1st Sort Filtering
        result = find_pattern_pair(sel_seq, s1_read1, s1_read2, max_mismatches, required_between=s1_length)
        if result is None:
            continue
        counts['sort1'] += 1

        # Stage 4: 2nd Sort Filtering
        result = find_pattern_pair(sel_seq, s2_read1, s2_read2, max_mismatches, required_between=s2_length)
        if result is None:
            continue
        counts['sort2'] += 1

        seq_counter[sel_seq] += 1

    return seq_counter, counts

def main():
    import argparse
    parser = argparse.ArgumentParser(description='AptaSelect: Aptamer candidate identification from SELEX data')
    parser.add_argument('--read1', required=True, help='Path to Read 1 FASTQ file')
    parser.add_argument('--read2', required=True, help='Path to Read 2 FASTQ file')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    parser.add_argument('--stats', required=True, help='Path to output stats JSON file')
    parser.add_argument('--config', required=True, help='Path to config JSON file')
    parser.add_argument('--chunk-size', type=int, default=10000, help='Chunk size for processing')
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    # Read paired FASTQ files
    total_counter = Counter()
    total_counts = {
        'total': 0,
        'joined': 0,
        'selection': 0,
        'sort1': 0,
        'sort2': 0,
    }

    chunk = []
    r1_iter = parse_fastq(args.read1)
    r2_iter = parse_fastq(args.read2)

    for (h1, s1, q1), (h2, s2, q2) in zip(r1_iter, r2_iter):
        chunk.append((s1, q1, s2, q2))
        if len(chunk) >= args.chunk_size:
            seq_counter, counts = process_chunk(chunk, config)
            total_counter += seq_counter
            for k in total_counts:
                total_counts[k] += counts[k]
            chunk = []

    # Process remaining records
    if chunk:
        seq_counter, counts = process_chunk(chunk, config)
        total_counter += seq_counter
        for k in total_counts:
            total_counts[k] += counts[k]

    # Stage 5: Aggregation & Ranking
    top_n = config.get('top_n', 10)
    all_sequences = total_counter.most_common()  # All unique sequences, sorted descending by count

    # Write output TSV — ALL unique sequences that passed the final filter
    with open(args.output, 'w') as f:
        f.write("Rank\tSequence\tCount\tLength\n")
        for rank, (seq, count) in enumerate(all_sequences, 1):
            f.write(f"{rank}\t{seq}\t{count}\t{len(seq)}\n")

    # Write stats JSON — includes top N summary for quick reference
    top_sequences = all_sequences[:top_n]
    stats = {
        'total_reads': total_counts['total'],
        'joined_reads': total_counts['joined'],
        'selection_passed': total_counts['selection'],
        'sort1_passed': total_counts['sort1'],
        'sort2_passed': total_counts['sort2'],
        'unique_sequences': len(total_counter),
        'top_n': top_n,
        'top_sequences': [{'rank': i+1, 'sequence': seq, 'count': cnt, 'length': len(seq)}
                          for i, (seq, cnt) in enumerate(top_sequences)]
    }
    with open(args.stats, 'w') as f:
        json.dump(stats, f, indent=2)

    # Print summary
    print(f"Total read pairs: {total_counts['total']}")
    print(f"Joined: {total_counts['joined']} ({100*total_counts['joined']/max(total_counts['total'],1):.1f}%)")
    print(f"Selection passed: {total_counts['selection']} ({100*total_counts['selection']/max(total_counts['total'],1):.1f}%)")
    print(f"1st Sort passed: {total_counts['sort1']} ({100*total_counts['sort1']/max(total_counts['total'],1):.1f}%)")
    print(f"2nd Sort passed: {total_counts['sort2']} ({100*total_counts['sort2']/max(total_counts['total'],1):.1f}%)")
    print(f"Unique sequences: {len(total_counter)}")
    print(f"All {len(all_sequences)} unique sequences written to {args.output}")

if __name__ == '__main__':
    main()