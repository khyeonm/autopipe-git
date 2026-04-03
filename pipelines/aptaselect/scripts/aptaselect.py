#!/usr/bin/env python3
"""
AptaSelect: Identify high-frequency aptamer candidate sequences from
paired-end FASTQ files produced by SELEX experiments.

Stages:
  1. Join paired-end reads via overlap detection
  2. Selection Filtering (extract region between left/right patterns)
  3. 1st Sort Filtering (validate patterns exist)
  4. 2nd Sort Filtering (validate patterns with between-length constraint)
  5. Aggregation & Ranking
"""

import argparse
import gzip
import os
import sys
import json
from collections import Counter
from itertools import islice


# ---------------------------------------------------------------------------
# FASTQ parsing
# ---------------------------------------------------------------------------

def parse_fastq(filepath):
    """Yield (name, seq, qual) tuples from a FASTQ file (plain or gzipped)."""
    opener = gzip.open if filepath.endswith(".gz") else open
    with opener(filepath, "rt") as fh:
        while True:
            header = fh.readline().rstrip("\n")
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()  # +
            qual = fh.readline().rstrip("\n")
            yield header, seq, qual


def reverse_complement(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def phred_score(ch):
    return ord(ch) - 33


# ---------------------------------------------------------------------------
# Stage 1: Join
# ---------------------------------------------------------------------------

def hamming(a, b):
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def join_reads(seq1, qual1, seq2, qual2, min_overlap=6, max_mismatch_pct=8.0,
               short_mode=True):
    """
    Join two reads by overlap.  Returns joined sequence or None.

    In short mode (default): first = rc(R2), second = R1
    In long mode:            first = R1,     second = rc(R2)
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
    best_score = None
    best_overlap = None

    for ov in range(min_overlap, max_overlap + 1):
        tail = first_seq[-ov:]
        head = second_seq[:ov]
        d = hamming(tail, head)
        rate = d / ov * 100.0
        if rate <= max_mismatch_pct:
            score = (1000 * (d * d + 1)) / ov
            if best_score is None or score < best_score:
                best_score = score
                best_overlap = ov

    if best_overlap is None:
        return None

    ov = best_overlap
    # Build merged overlap region using higher-quality base
    overlap_seq = []
    tail_start = len(first_seq) - ov
    for i in range(ov):
        q1 = phred_score(first_qual[tail_start + i])
        q2 = phred_score(second_qual[i])
        if q1 >= q2:
            overlap_seq.append(first_seq[tail_start + i])
        else:
            overlap_seq.append(second_seq[i])

    joined = first_seq[:tail_start] + "".join(overlap_seq) + second_seq[ov:]
    return joined


# ---------------------------------------------------------------------------
# Stages 2-4: Pattern matching helpers
# ---------------------------------------------------------------------------

def find_pattern(seq, pattern, max_mismatch=1):
    """Yield all start positions where pattern matches seq within tolerance."""
    plen = len(pattern)
    for i in range(len(seq) - plen + 1):
        d = hamming(seq[i:i + plen], pattern)
        if d <= max_mismatch:
            yield i


def find_pattern_pair(seq, left_pat, right_pat, max_mismatch=1,
                      between_length=None):
    """
    Nested-loop search for left then right pattern.
    Returns (left_start, right_start) or None.
    """
    lplen = len(left_pat)
    rplen = len(right_pat)
    for lpos in find_pattern(seq, left_pat, max_mismatch):
        region_start = lpos + lplen  # after left pattern
        for rpos in find_pattern(seq[region_start:], right_pat, max_mismatch):
            abs_rpos = region_start + rpos
            between = abs_rpos - region_start
            if between_length is not None and between != between_length:
                continue
            return lpos, abs_rpos
    return None


# ---------------------------------------------------------------------------
# Stage 2: Selection Filtering
# ---------------------------------------------------------------------------

def selection_filter(seq, left_pat, right_pat, max_mismatch=1):
    """Extract region from start-of-left to end-of-right (inclusive)."""
    result = find_pattern_pair(seq, left_pat, right_pat, max_mismatch)
    if result is None:
        return None
    lpos, rpos = result
    return seq[lpos: rpos + len(right_pat)]


# ---------------------------------------------------------------------------
# Stage 3: 1st Sort Filtering
# ---------------------------------------------------------------------------

def sort1_filter(seq, left_pat, right_pat, max_mismatch=1):
    """Return seq unchanged if both patterns found, else None."""
    result = find_pattern_pair(seq, left_pat, right_pat, max_mismatch)
    if result is None:
        return None
    return seq


# ---------------------------------------------------------------------------
# Stage 4: 2nd Sort Filtering
# ---------------------------------------------------------------------------

def sort2_filter(seq, left_pat, right_pat, max_mismatch=1,
                 between_length=20):
    """Return seq unchanged if patterns found with exact between-length."""
    result = find_pattern_pair(seq, left_pat, right_pat, max_mismatch,
                               between_length=between_length)
    if result is None:
        return None
    return seq


# ---------------------------------------------------------------------------
# Stage 5: Aggregation helpers
# ---------------------------------------------------------------------------

def write_ranked(counter, filepath):
    """Write frequency-ranked sequences to a TSV file."""
    with open(filepath, "w") as fh:
        fh.write("rank\tcount\tsequence\n")
        for rank, (seq, cnt) in enumerate(
                counter.most_common(), start=1):
            fh.write(f"{rank}\t{cnt}\t{seq}\n")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def chunked(iterable, n):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, n))
        if not chunk:
            return
        yield chunk


def run_pipeline(args):
    # Counters for each stage
    join_counter = Counter()
    sel_counter = Counter()
    sort1_counter = Counter()
    sort2_counter = Counter()

    stats = {
        "total_pairs": 0,
        "joined": 0,
        "selection_pass": 0,
        "sort1_pass": 0,
        "sort2_pass": 0,
    }

    r1_iter = parse_fastq(args.read1)
    r2_iter = parse_fastq(args.read2)

    chunk_idx = 0
    for chunk in chunked(zip(r1_iter, r2_iter), args.chunk_size):
        chunk_idx += 1
        for (h1, s1, q1), (h2, s2, q2) in chunk:
            stats["total_pairs"] += 1

            # Stage 1: Join
            joined = join_reads(
                s1, q1, s2, q2,
                min_overlap=args.min_overlap,
                max_mismatch_pct=args.max_mismatch_pct,
                short_mode=args.short_mode,
            )
            if joined is None:
                continue
            stats["joined"] += 1
            join_counter[joined] += 1

            # Stage 2: Selection Filtering
            extracted = selection_filter(
                joined,
                args.sel_left, args.sel_right,
                max_mismatch=args.pattern_mismatch,
            )
            if extracted is None:
                continue
            stats["selection_pass"] += 1
            sel_counter[extracted] += 1

            # Stage 3: 1st Sort Filtering
            s3 = sort1_filter(
                extracted,
                args.sort1_left, args.sort1_right,
                max_mismatch=args.pattern_mismatch,
            )
            if s3 is None:
                continue
            stats["sort1_pass"] += 1
            sort1_counter[s3] += 1

            # Stage 4: 2nd Sort Filtering
            s4 = sort2_filter(
                extracted,
                args.sort2_left, args.sort2_right,
                max_mismatch=args.pattern_mismatch,
                between_length=args.sort2_between_length,
            )
            if s4 is None:
                continue
            stats["sort2_pass"] += 1
            sort2_counter[s4] += 1

    # Stage 5: Write results
    os.makedirs(args.outdir, exist_ok=True)

    write_ranked(join_counter, os.path.join(args.outdir, "stage1_joined.tsv"))
    write_ranked(sel_counter, os.path.join(args.outdir, "stage2_selection.tsv"))
    write_ranked(sort1_counter, os.path.join(args.outdir, "stage3_sort1.tsv"))
    write_ranked(sort2_counter, os.path.join(args.outdir, "stage4_sort2.tsv"))

    # Write stats
    stats_path = os.path.join(args.outdir, "pipeline_stats.json")
    with open(stats_path, "w") as fh:
        json.dump(stats, fh, indent=2)

    print("=== AptaSelect Pipeline Complete ===")
    print(f"Total read pairs:       {stats['total_pairs']}")
    print(f"Joined (Stage 1):       {stats['joined']}")
    print(f"Selection pass (St 2):  {stats['selection_pass']}")
    print(f"Sort-1 pass (St 3):     {stats['sort1_pass']}")
    print(f"Sort-2 pass (St 4):     {stats['sort2_pass']}")
    print(f"Unique seqs (Stage 2):  {len(sel_counter)}")
    print(f"Unique seqs (Stage 3):  {len(sort1_counter)}")
    print(f"Unique seqs (Stage 4):  {len(sort2_counter)}")


def main():
    p = argparse.ArgumentParser(description="AptaSelect aptamer pipeline")

    # I/O
    p.add_argument("--read1", required=True, help="Path to Read 1 FASTQ")
    p.add_argument("--read2", required=True, help="Path to Read 2 FASTQ")
    p.add_argument("--outdir", required=True, help="Output directory")

    # Stage 1 params
    p.add_argument("--min-overlap", type=int, default=6)
    p.add_argument("--max-mismatch-pct", type=float, default=8.0)
    p.add_argument("--short-mode", action="store_true", default=True,
                   help="Short library mode (default: enabled)")
    p.add_argument("--long-mode", dest="short_mode", action="store_false",
                   help="Long library mode")

    # Pattern mismatch tolerance (shared across stages 2-4)
    p.add_argument("--pattern-mismatch", type=int, default=1)

    # Stage 2 patterns
    p.add_argument("--sel-left", default="CCACTTCTCCTTCCATCCTAAAC")
    p.add_argument("--sel-right", default="GAGTAGTTTGGAGGGTTGTCTG")

    # Stage 3 patterns
    p.add_argument("--sort1-left", default="TCCTAAAC")
    p.add_argument("--sort1-right", default="GAGTAGTT")

    # Stage 4 patterns
    p.add_argument("--sort2-left", default="TCTCTCTCTC")
    p.add_argument("--sort2-right", default="GAGAGAGAGA")
    p.add_argument("--sort2-between-length", type=int, default=20)

    # Processing
    p.add_argument("--chunk-size", type=int, default=10000)

    args = p.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()