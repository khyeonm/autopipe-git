#!/usr/bin/env python3
"""
AptaSelect: Identifies high-frequency aptamer candidate sequences from paired-end FASTQ files
produced by SELEX experiments.

Stages:
  1. Join - Overlap-merge paired-end reads
  2. Selection Filtering - Extract region between selection patterns
  3. 1st Sort Filtering - Validate presence of 1st sort patterns
  4. 2nd Sort Filtering - Validate presence of 2nd sort patterns with length constraint
  5. Aggregation & Ranking - Deduplicate, count, and rank sequences
"""

import argparse
import gzip
import sys
import os
from collections import defaultdict, Counter


# ──────────────────────────────────────────────────────────────────────────────
# FASTQ utilities
# ──────────────────────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fastq_pairs(r1_path, r2_path):
    """Yield (name, seq1, qual1, seq2, qual2) for each read pair."""
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2:
        while True:
            # Read 1
            name1 = f1.readline().rstrip()
            if not name1:
                break
            seq1 = f1.readline().rstrip()
            f1.readline()  # +
            qual1 = f1.readline().rstrip()
            # Read 2
            name2 = f2.readline().rstrip()
            seq2 = f2.readline().rstrip()
            f2.readline()  # +
            qual2 = f2.readline().rstrip()
            yield (name1, seq1, qual1, seq2, qual2)


# ──────────────────────────────────────────────────────────────────────────────
# Stage 1: Join
# ──────────────────────────────────────────────────────────────────────────────

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def join_reads(seq1, qual1, seq2, qual2, min_overlap=6, max_mismatch_pct=0.08,
               short_mode=True):
    """
    Join paired-end reads by finding optimal overlap.
    Returns joined sequence or None if no valid overlap found.
    """
    rc_seq2 = reverse_complement(seq2)
    rc_qual2 = qual2[::-1]

    if short_mode:
        first_seq, first_qual = rc_seq2, rc_qual2
        second_seq, second_qual = seq1, qual1
    else:
        first_seq, first_qual = seq1, qual1
        second_seq, second_qual = rc_seq2, rc_qual2

    max_overlap = min(len(first_seq), len(second_seq))
    best_score = float('inf')
    best_overlap = -1

    for ov in range(min_overlap, max_overlap + 1):
        tail = first_seq[len(first_seq) - ov:]
        head = second_seq[:ov]
        d = hamming_distance(tail, head)
        mismatch_rate = d / ov
        if mismatch_rate <= max_mismatch_pct:
            score = (1000 * (d * d + 1)) / ov
            if score < best_score:
                best_score = score
                best_overlap = ov

    if best_overlap < 0:
        return None

    # Resolve mismatches using quality scores
    ov = best_overlap
    first_tail_start = len(first_seq) - ov
    overlap_seq = []
    for i in range(ov):
        fi = first_tail_start + i
        si = i
        if first_seq[fi] == second_seq[si]:
            overlap_seq.append(first_seq[fi])
        else:
            q1 = ord(first_qual[fi]) - 33
            q2 = ord(second_qual[si]) - 33
            overlap_seq.append(first_seq[fi] if q1 >= q2 else second_seq[si])

    joined = first_seq[:first_tail_start] + "".join(overlap_seq) + second_seq[ov:]
    return joined


# ──────────────────────────────────────────────────────────────────────────────
# Stages 2–4: Pattern matching
# ──────────────────────────────────────────────────────────────────────────────

def pattern_matches(seq, pattern, max_mismatches):
    """Yield all start positions where pattern matches seq within mismatch tolerance."""
    plen = len(pattern)
    for i in range(len(seq) - plen + 1):
        d = hamming_distance(seq[i:i + plen], pattern)
        if d <= max_mismatches:
            yield i


def find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches=1,
                      required_between_length=None):
    """
    Nested-loop search for left and right pattern positions.
    Returns (left_start, right_end) or None.
    """
    lp_len = len(left_pattern)
    rp_len = len(right_pattern)

    for left_pos in pattern_matches(seq, left_pattern, max_mismatches):
        region_start = left_pos + lp_len
        for right_pos in pattern_matches(seq[region_start:], right_pattern, max_mismatches):
            abs_right_pos = region_start + right_pos
            between_len = abs_right_pos - region_start
            if required_between_length is not None and between_len != required_between_length:
                continue
            right_end = abs_right_pos + rp_len
            return (left_pos, right_end)
    return None


def stage2_selection(seq, left_pattern, right_pattern, max_mismatches):
    """Extract the region from start of left pattern to end of right pattern."""
    result = find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches)
    if result is None:
        return None
    left_start, right_end = result
    return seq[left_start:right_end]


def stage3_sort1(seq, left_pattern, right_pattern, max_mismatches):
    """Validate 1st sort patterns exist; pass through unchanged if valid."""
    result = find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches)
    return seq if result is not None else None


def stage4_sort2(seq, left_pattern, right_pattern, max_mismatches, between_length):
    """Validate 2nd sort patterns with required between-length; pass through unchanged."""
    result = find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches,
                               required_between_length=between_length)
    return seq if result is not None else None


# ──────────────────────────────────────────────────────────────────────────────
# Stage 5: Aggregation & Ranking
# ──────────────────────────────────────────────────────────────────────────────

def aggregate_and_rank(counter):
    """Return list of (sequence, count) sorted by count descending."""
    return sorted(counter.items(), key=lambda x: -x[1])


# ──────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ──────────────────────────────────────────────────────────────────────────────

def write_ranked(filepath, ranked):
    with open(filepath, "w") as f:
        f.write("rank\tcount\tsequence\n")
        for rank, (seq, count) in enumerate(ranked, 1):
            f.write(f"{rank}\t{count}\t{seq}\n")


def run_pipeline(args):
    # Counters for each stage
    stage1_counter = Counter()
    stage2_counter = Counter()
    stage3_counter = Counter()
    stage4_counter = Counter()

    total_pairs = 0
    joined_count = 0
    s2_pass = 0
    s3_pass = 0
    s4_pass = 0

    chunk_size = args.chunk_size
    chunk_seqs = []

    def process_chunk(seqs):
        nonlocal joined_count, s2_pass, s3_pass, s4_pass
        for seq in seqs:
            # Stage 1 is already done; seq is the joined sequence
            stage1_counter[seq] += 1

            # Stage 2: Selection Filtering
            extracted = stage2_selection(seq, args.sel_left, args.sel_right,
                                         args.max_mismatches)
            if extracted is None:
                continue
            s2_pass += 1
            stage2_counter[extracted] += 1

            # Stage 3: 1st Sort Filtering
            passed = stage3_sort1(extracted, args.sort1_left, args.sort1_right,
                                  args.max_mismatches)
            if passed is None:
                continue
            s3_pass += 1
            stage3_counter[passed] += 1

            # Stage 4: 2nd Sort Filtering
            passed = stage4_sort2(extracted, args.sort2_left, args.sort2_right,
                                  args.max_mismatches, args.sort2_between_length)
            if passed is None:
                continue
            s4_pass += 1
            stage4_counter[passed] += 1

    print(f"Processing {args.r1} and {args.r2} ...")
    print(f"Short mode: {args.short_mode}")

    for name, seq1, qual1, seq2, qual2 in read_fastq_pairs(args.r1, args.r2):
        total_pairs += 1

        # Stage 1: Join
        joined = join_reads(seq1, qual1, seq2, qual2,
                           min_overlap=args.min_overlap,
                           max_mismatch_pct=args.max_mismatch_pct,
                           short_mode=args.short_mode)
        if joined is None:
            continue
        joined_count += 1
        chunk_seqs.append(joined)

        if len(chunk_seqs) >= chunk_size:
            process_chunk(chunk_seqs)
            chunk_seqs = []

    # Process remaining
    if chunk_seqs:
        process_chunk(chunk_seqs)

    # Summary
    print(f"\n=== AptaSelect Summary ===")
    print(f"Total read pairs:           {total_pairs}")
    print(f"Stage 1 (Join):             {joined_count} ({joined_count/max(total_pairs,1)*100:.1f}%)")
    print(f"Stage 2 (Selection Filter): {s2_pass} ({s2_pass/max(total_pairs,1)*100:.1f}%)")
    print(f"Stage 3 (1st Sort Filter):  {s3_pass} ({s3_pass/max(total_pairs,1)*100:.1f}%)")
    print(f"Stage 4 (2nd Sort Filter):  {s4_pass} ({s4_pass/max(total_pairs,1)*100:.1f}%)")

    # Stage 5: Aggregate & Rank each stage
    os.makedirs(args.outdir, exist_ok=True)

    ranked1 = aggregate_and_rank(stage1_counter)
    write_ranked(os.path.join(args.outdir, "stage1_joined_ranked.tsv"), ranked1)
    print(f"Stage 1 unique sequences:   {len(ranked1)}")

    ranked2 = aggregate_and_rank(stage2_counter)
    write_ranked(os.path.join(args.outdir, "stage2_selection_ranked.tsv"), ranked2)
    print(f"Stage 2 unique sequences:   {len(ranked2)}")

    ranked3 = aggregate_and_rank(stage3_counter)
    write_ranked(os.path.join(args.outdir, "stage3_sort1_ranked.tsv"), ranked3)
    print(f"Stage 3 unique sequences:   {len(ranked3)}")

    ranked4 = aggregate_and_rank(stage4_counter)
    write_ranked(os.path.join(args.outdir, "stage4_sort2_ranked.tsv"), ranked4)
    print(f"Stage 4 unique sequences:   {len(ranked4)}")

    # Write summary
    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write(f"Total read pairs:           {total_pairs}\n")
        f.write(f"Stage 1 (Join):             {joined_count} ({joined_count/max(total_pairs,1)*100:.1f}%)\n")
        f.write(f"Stage 2 (Selection Filter): {s2_pass} ({s2_pass/max(total_pairs,1)*100:.1f}%)\n")
        f.write(f"Stage 3 (1st Sort Filter):  {s3_pass} ({s3_pass/max(total_pairs,1)*100:.1f}%)\n")
        f.write(f"Stage 4 (2nd Sort Filter):  {s4_pass} ({s4_pass/max(total_pairs,1)*100:.1f}%)\n")
        f.write(f"Stage 1 unique sequences:   {len(ranked1)}\n")
        f.write(f"Stage 2 unique sequences:   {len(ranked2)}\n")
        f.write(f"Stage 3 unique sequences:   {len(ranked3)}\n")
        f.write(f"Stage 4 unique sequences:   {len(ranked4)}\n")

    print("\nDone! Results written to:", args.outdir)


def main():
    parser = argparse.ArgumentParser(description="AptaSelect: SELEX aptamer candidate identification")

    # Input
    parser.add_argument("--r1", required=True, help="Read 1 FASTQ file path")
    parser.add_argument("--r2", required=True, help="Read 2 FASTQ file path")
    parser.add_argument("--outdir", required=True, help="Output directory")

    # Stage 1: Join parameters
    parser.add_argument("--short-mode", action="store_true", default=True,
                        help="Use short library mode (default: True)")
    parser.add_argument("--long-mode", action="store_true", default=False,
                        help="Use long library mode (overrides short-mode)")
    parser.add_argument("--min-overlap", type=int, default=6,
                        help="Minimum overlap length in bp (default: 6)")
    parser.add_argument("--max-mismatch-pct", type=float, default=0.08,
                        help="Maximum mismatch rate for overlap (default: 0.08)")

    # Pattern matching tolerance
    parser.add_argument("--max-mismatches", type=int, default=1,
                        help="Maximum mismatches for pattern matching (default: 1)")

    # Stage 2: Selection patterns
    parser.add_argument("--sel-left", default="CCACTTCTCCTTCCATCCTAAAC",
                        help="Selection left pattern")
    parser.add_argument("--sel-right", default="GAGTAGTTTGGAGGGTTGTCTG",
                        help="Selection right pattern")

    # Stage 3: 1st Sort patterns
    parser.add_argument("--sort1-left", default="TCCTAAAC",
                        help="1st Sort left pattern")
    parser.add_argument("--sort1-right", default="GAGTAGTT",
                        help="1st Sort right pattern")

    # Stage 4: 2nd Sort patterns
    parser.add_argument("--sort2-left", default="TCTCTCTCTC",
                        help="2nd Sort left pattern")
    parser.add_argument("--sort2-right", default="GAGAGAGAGA",
                        help="2nd Sort right pattern")
    parser.add_argument("--sort2-between-length", type=int, default=20,
                        help="Required between-length for 2nd Sort (default: 20)")

    # Chunk size
    parser.add_argument("--chunk-size", type=int, default=10000,
                        help="Processing chunk size (default: 10000)")

    args = parser.parse_args()

    if args.long_mode:
        args.short_mode = False

    run_pipeline(args)


if __name__ == "__main__":
    main()