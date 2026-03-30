#!/usr/bin/env python3
"""
AptaSelect: Identify high-frequency aptamer candidate sequences from paired-end
FASTQ files produced by SELEX experiments.

Pipeline stages:
  1. Join paired-end reads via overlap detection
  2. Selection Filtering (extract region between selection patterns)
  3. 1st Sort Filtering (validate between-length for 1st sort patterns)
  4. 2nd Sort Filtering (validate between-length for 2nd sort patterns)
  5. Aggregation & Ranking (count, sort, report all unique sequences)
"""

import argparse
import gzip
import json
import sys
from collections import Counter


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def phred_char_to_score(ch):
    """Convert a Phred+33 ASCII character to an integer quality score."""
    return ord(ch) - 33


def parse_fastq(path):
    """Yield (name, seq, qual) tuples from a FASTQ file (plain or gzipped)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline().rstrip("\n")
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # + line
            qual = fh.readline().rstrip("\n")
            yield header, seq, qual


# ---------------------------------------------------------------------------
# Stage 1: Join paired-end reads
# ---------------------------------------------------------------------------

def hamming(a, b):
    """Return the Hamming distance between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def join_reads(seq1, qual1, seq2, qual2, min_overlap=6, max_mismatch_rate=0.08,
               short_mode=True):
    """
    Join two reads by finding the best overlap.

    In short mode (default): RC(R2) is the first sequence, R1 is the second.
    In long mode: R1 is the first sequence, RC(R2) is the second.

    Returns the joined sequence or None if no valid overlap is found.
    """
    rc_seq2 = reverse_complement(seq2)
    rc_qual2 = qual2[::-1]

    if short_mode:
        first_seq, first_qual = rc_seq2, rc_qual2
        second_seq, second_qual = seq1, qual1
    else:
        first_seq, first_qual = seq1, qual1
        second_seq, second_qual = rc_seq2, rc_qual2

    best_score = float("inf")
    best_overlap = None

    max_overlap = min(len(first_seq), len(second_seq))

    for ovl in range(min_overlap, max_overlap + 1):
        tail = first_seq[-ovl:]
        head = second_seq[:ovl]
        d = hamming(tail, head)
        rate = d / ovl
        if rate <= max_mismatch_rate:
            score = (1000 * (d * d + 1)) / ovl
            if score < best_score:
                best_score = score
                best_overlap = ovl

    if best_overlap is None:
        return None

    ovl = best_overlap
    # Build the overlap consensus using quality scores
    overlap_seq = []
    for i in range(ovl):
        pos_first = len(first_seq) - ovl + i
        pos_second = i
        if first_seq[pos_first] == second_seq[pos_second]:
            overlap_seq.append(first_seq[pos_first])
        else:
            q1 = phred_char_to_score(first_qual[pos_first])
            q2 = phred_char_to_score(second_qual[pos_second])
            overlap_seq.append(first_seq[pos_first] if q1 >= q2 else second_seq[pos_second])

    joined = first_seq[:-ovl] + "".join(overlap_seq) + second_seq[ovl:]
    return joined


# ---------------------------------------------------------------------------
# Stages 2-4: Pattern-based filtering
# ---------------------------------------------------------------------------

def find_pattern(seq, pattern, max_mismatches=1):
    """Yield all start positions where pattern matches seq within tolerance."""
    plen = len(pattern)
    for i in range(len(seq) - plen + 1):
        if hamming(seq[i:i + plen], pattern) <= max_mismatches:
            yield i


def pattern_filter(seq, left_pattern, right_pattern, max_mismatches=1,
                   required_between_length=None, extract=False):
    """
    Search for left and right patterns using nested loops.

    If extract=True, return the region from start-of-left to end-of-right.
    Otherwise, validate and return the input sequence unchanged.

    Returns the (possibly extracted) sequence, or None if no valid pair found.
    """
    left_len = len(left_pattern)
    right_len = len(right_pattern)

    for left_pos in find_pattern(seq, left_pattern, max_mismatches):
        region_start = left_pos + left_len  # position right after left pattern
        for right_pos in find_pattern(seq[region_start:], right_pattern, max_mismatches):
            abs_right_pos = region_start + right_pos  # absolute position of right pattern start
            between_len = abs_right_pos - region_start
            if required_between_length is not None and between_len != required_between_length:
                continue
            # Valid pair found
            if extract:
                return seq[left_pos: abs_right_pos + right_len]
            else:
                return seq
    return None


# ---------------------------------------------------------------------------
# Stage 5: Aggregation & Ranking
# ---------------------------------------------------------------------------

def aggregate_and_rank(sequences, chunk_size=10000):
    """Count unique sequences in chunks, merge, and return all sorted by count descending."""
    total_counter = Counter()
    chunk = []

    for s in sequences:
        chunk.append(s)
        if len(chunk) >= chunk_size:
            total_counter.update(chunk)
            chunk = []

    if chunk:
        total_counter.update(chunk)

    return total_counter.most_common()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(r1_path, r2_path, output_path, stats_path, cfg):
    """Execute the full AptaSelect pipeline."""

    # Unpack config
    min_overlap = cfg.get("min_overlap", 6)
    max_mismatch_rate = cfg.get("max_mismatch_rate", 0.08)
    short_mode = cfg.get("short_mode", True)
    pattern_mismatches = cfg.get("pattern_mismatches", 1)

    sel_left = cfg.get("selection_left_pattern", "CCACTTCTCCTTCCATCCTAAAC")
    sel_right = cfg.get("selection_right_pattern", "GAGTAGTTTGGAGGGTTGTCTG")
    sort1_left = cfg.get("sort1_left_pattern", "TCCTAAAC")
    sort1_right = cfg.get("sort1_right_pattern", "GAGTAGTT")
    sort1_between = cfg.get("sort1_between_length", 40)
    sort2_left = cfg.get("sort2_left_pattern", "TCTCTCTCTC")
    sort2_right = cfg.get("sort2_right_pattern", "GAGAGAGAGA")
    sort2_between = cfg.get("sort2_between_length", 20)

    chunk_size = cfg.get("chunk_size", 10000)

    # Counters for stats
    total_pairs = 0
    joined_count = 0
    stage2_count = 0
    stage3_count = 0
    stage4_count = 0

    passed_sequences = []

    for (h1, s1, q1), (h2, s2, q2) in zip(parse_fastq(r1_path), parse_fastq(r2_path)):
        total_pairs += 1

        # Stage 1: Join
        joined = join_reads(s1, q1, s2, q2,
                            min_overlap=min_overlap,
                            max_mismatch_rate=max_mismatch_rate,
                            short_mode=short_mode)
        if joined is None:
            continue
        joined_count += 1

        # Stage 2: Selection Filtering (extract)
        extracted = pattern_filter(joined, sel_left, sel_right,
                                   max_mismatches=pattern_mismatches,
                                   required_between_length=None,
                                   extract=True)
        if extracted is None:
            continue
        stage2_count += 1

        # Stage 3: 1st Sort Filtering (validate, pass-through)
        result3 = pattern_filter(extracted, sort1_left, sort1_right,
                                  max_mismatches=pattern_mismatches,
                                  required_between_length=sort1_between,
                                  extract=False)
        if result3 is None:
            continue
        stage3_count += 1

        # Stage 4: 2nd Sort Filtering (validate, pass-through)
        result4 = pattern_filter(extracted, sort2_left, sort2_right,
                                  max_mismatches=pattern_mismatches,
                                  required_between_length=sort2_between,
                                  extract=False)
        if result4 is None:
            continue
        stage4_count += 1

        passed_sequences.append(extracted)

    # Stage 5: Aggregation & Ranking (all unique sequences)
    ranked = aggregate_and_rank(passed_sequences, chunk_size=chunk_size)

    # Write ranked output (all unique sequences, sorted by count descending)
    with open(output_path, "w") as fh:
        fh.write("Rank\tSequence\tCount\n")
        for rank, (seq, count) in enumerate(ranked, 1):
            fh.write(f"{rank}\t{seq}\t{count}\n")

    # Write stats
    stats = {
        "total_read_pairs": total_pairs,
        "stage1_joined": joined_count,
        "stage2_selection_passed": stage2_count,
        "stage3_sort1_passed": stage3_count,
        "stage4_sort2_passed": stage4_count,
        "stage5_unique_sequences": len(ranked),
        "stage5_total_passed": len(passed_sequences),
    }
    with open(stats_path, "w") as fh:
        json.dump(stats, fh, indent=2)

    # Print summary
    print(f"Total read pairs:        {total_pairs}")
    print(f"Stage 1 (Joined):        {joined_count}")
    print(f"Stage 2 (Selection):     {stage2_count}")
    print(f"Stage 3 (1st Sort):      {stage3_count}")
    print(f"Stage 4 (2nd Sort):      {stage4_count}")
    print(f"Unique sequences:        {len(ranked)}")
    print(f"Total passed:            {len(passed_sequences)}")
    print(f"All {len(ranked)} unique candidates written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AptaSelect pipeline")
    parser.add_argument("--r1", required=True, help="Read 1 FASTQ file")
    parser.add_argument("--r2", required=True, help="Read 2 FASTQ file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--stats", required=True, help="Stats JSON file")
    parser.add_argument("--config", required=True, help="Config YAML file")
    args = parser.parse_args()

    import yaml
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    run_pipeline(args.r1, args.r2, args.output, args.stats, cfg)