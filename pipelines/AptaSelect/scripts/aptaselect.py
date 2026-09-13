#!/usr/bin/env python3
"""
AptaSelect: identify high-frequency aptamer candidate sequences from paired-end
FASTQ files produced by SELEX experiments.

Pipeline stages (each read mate processed independently, survivors pooled):
  Stage 0 - Raw: aggregate all raw reads by identical sequence.
  Stage 1 - Matching: orient each read by LCR...RCR (mismatch-tolerant), checking
            both the read and its reverse complement.
  Stage 2 - Flanking extraction: extract the region between an inner LCR marker
            and an inner RCR marker (exact match, zero mismatch).
  Stage 3 - N extraction: within the extracted region, find left/right stems by
            exact match and require the between-region to be exactly the target
            random-region length.

Every pattern sequence, mismatch tolerance and length is passed in from
config.yaml (nothing hard-coded here).
"""

import argparse
import gzip
import sys
from collections import OrderedDict


# --------------------------------------------------------------------------- #
# FASTQ reading
# --------------------------------------------------------------------------- #
def open_maybe_gzip(path):
    """Open a file transparently whether it is gzip-compressed or plain text."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fastq_sequences(path):
    """Yield the sequence line (2nd line of every 4-line record) from a FASTQ."""
    with open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break  # truncated record at EOF
            seq = seq.strip()
            if seq:
                yield seq


# --------------------------------------------------------------------------- #
# Sequence utilities
# --------------------------------------------------------------------------- #
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def hamming_at_most(text, pattern, start, max_mismatch):
    """True if pattern matches text starting at `start` within max_mismatch."""
    mm = 0
    for i in range(len(pattern)):
        if text[start + i] != pattern[i]:
            mm += 1
            if mm > max_mismatch:
                return False
    return True


def find_first_fuzzy(text, pattern, max_mismatch, from_pos=0):
    """
    Return the start index of the first occurrence of `pattern` in `text`
    (scanning left to right, at or after from_pos) that matches within
    max_mismatch substitutions, else -1.
    """
    plen = len(pattern)
    tlen = len(text)
    if plen == 0:
        return from_pos
    last_start = tlen - plen
    for start in range(from_pos, last_start + 1):
        if hamming_at_most(text, pattern, start, max_mismatch):
            return start
    return -1


def find_all_exact(text, pattern):
    """Return list of all start indices of exact (possibly overlapping) matches."""
    positions = []
    plen = len(pattern)
    if plen == 0:
        return positions
    start = text.find(pattern)
    while start != -1:
        positions.append(start)
        start = text.find(pattern, start + 1)
    return positions


# --------------------------------------------------------------------------- #
# Stage 1 - Matching (orient by LCR ... RCR, mismatch-tolerant, both strands)
# --------------------------------------------------------------------------- #
def orient_read(seq, lcr, rcr, max_mismatch):
    """
    Return the read oriented so that LCR lies upstream of RCR, or None if the
    LCR...RCR pair is found in neither the read nor its reverse complement.

    An LCR...RCR pair is "found" when an LCR match exists and an RCR match
    exists that begins at or after the end of that LCR match.
    """
    for candidate in (seq, reverse_complement(seq)):
        lcr_start = find_first_fuzzy(candidate, lcr, max_mismatch, 0)
        if lcr_start == -1:
            continue
        lcr_end = lcr_start + len(lcr)
        rcr_start = find_first_fuzzy(candidate, rcr, max_mismatch, lcr_end)
        if rcr_start != -1:
            return candidate
    return None


# --------------------------------------------------------------------------- #
# Stage 2 - Flanking extraction (between inner markers, exact match)
# --------------------------------------------------------------------------- #
def extract_flanking(seq, left_marker, right_marker):
    """
    Find left_marker (inner end of LCR) and right_marker (inner start of RCR)
    by exact match, and return the sequence lying strictly between them.
    Return None if either marker is not found.

    The left marker is taken as its first exact occurrence; the right marker is
    taken as its first exact occurrence that begins at or after the end of the
    left marker (so the extracted region is well defined and non-crossing).
    """
    lpos = seq.find(left_marker)
    if lpos == -1:
        return None
    inner_start = lpos + len(left_marker)
    rpos = seq.find(right_marker, inner_start)
    if rpos == -1:
        return None
    return seq[inner_start:rpos]


# --------------------------------------------------------------------------- #
# Stage 3 - N (random region) extraction (between stems, exact target length)
# --------------------------------------------------------------------------- #
def extract_random_region(seq, left_stem, right_stem, target_len):
    """
    Find left_stem and right_stem by exact match and return the region between
    them when it is exactly target_len long. When a repetitive stem matches at
    several overlapping positions, consider all left/right position pairs and
    accept the FIRST pair (left positions ascending, then right positions
    ascending) that satisfies the exact-length constraint. Return None if no
    such pair exists.
    """
    left_positions = find_all_exact(seq, left_stem)
    if not left_positions:
        return None
    right_positions = find_all_exact(seq, right_stem)
    if not right_positions:
        return None

    llen = len(left_stem)
    for lpos in left_positions:
        inner_start = lpos + llen
        for rpos in right_positions:
            if rpos < inner_start:
                continue
            region = seq[inner_start:rpos]
            if len(region) == target_len:
                return region
    return None


# --------------------------------------------------------------------------- #
# Aggregation and ranking
# --------------------------------------------------------------------------- #
def aggregate_and_write(sequences, out_path):
    """
    Aggregate identical sequences by count and write a ranked TSV
    (sequence<TAB>count), sorted by count descending with first-seen order as
    the stable tie-breaker. Returns (total, unique).
    """
    counts = OrderedDict()  # preserves first-seen order
    total = 0
    for s in sequences:
        total += 1
        counts[s] = counts.get(s, 0) + 1

    # Stable sort by -count keeps first-seen order among equal counts.
    items = sorted(counts.items(), key=lambda kv: -kv[1])

    with open(out_path, "w") as out:
        out.write("sequence\tcount\n")
        for seq, cnt in items:
            out.write(f"{seq}\t{cnt}\n")

    return total, len(counts)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="AptaSelect SELEX aptamer filter")
    ap.add_argument("--r1", required=True, help="Read 1 FASTQ (.fq/.fastq[.gz])")
    ap.add_argument("--r2", required=True, help="Read 2 FASTQ (.fq/.fastq[.gz])")
    ap.add_argument("--outdir", required=True, help="Output directory")

    ap.add_argument("--lcr", required=True, help="Left constant region pattern")
    ap.add_argument("--rcr", required=True, help="Right constant region pattern")
    ap.add_argument("--constant-mismatch", type=int, required=True,
                    help="Max mismatches allowed for LCR / RCR matching")

    ap.add_argument("--left-marker", required=True,
                    help="Inner-end-of-LCR marker (exact match)")
    ap.add_argument("--right-marker", required=True,
                    help="Inner-start-of-RCR marker (exact match)")

    ap.add_argument("--left-stem", required=True, help="Left stem (exact match)")
    ap.add_argument("--right-stem", required=True, help="Right stem (exact match)")
    ap.add_argument("--random-region-length", type=int, required=True,
                    help="Exact length of the random region N")

    args = ap.parse_args()

    outdir = args.outdir.rstrip("/")

    # Normalize sequence patterns to uppercase for consistent matching.
    lcr = args.lcr.strip().upper()
    rcr = args.rcr.strip().upper()
    left_marker = args.left_marker.strip().upper()
    right_marker = args.right_marker.strip().upper()
    left_stem = args.left_stem.strip().upper()
    right_stem = args.right_stem.strip().upper()
    const_mm = args.constant_mismatch
    n_len = args.random_region_length

    # ------------------------------------------------------------------ #
    # Read both mates. Each mate is its own sequence; survivors are pooled.
    # ------------------------------------------------------------------ #
    raw_reads = []
    for path in (args.r1, args.r2):
        for seq in read_fastq_sequences(path):
            raw_reads.append(seq.upper())

    # Stage 0 - raw
    total_reads, raw_unique = aggregate_and_write(
        raw_reads, f"{outdir}/stage0_raw.ranked.tsv"
    )

    # Stage 1 - matching / orientation
    stage1 = []
    for seq in raw_reads:
        oriented = orient_read(seq, lcr, rcr, const_mm)
        if oriented is not None:
            stage1.append(oriented)
    aggregate_and_write(stage1, f"{outdir}/stage1_matching.ranked.tsv")

    # Stage 2 - flanking extraction (carries forward the oriented sequence)
    stage2 = []
    for seq in stage1:
        region = extract_flanking(seq, left_marker, right_marker)
        if region is not None:
            stage2.append(region)
    aggregate_and_write(stage2, f"{outdir}/stage2_flanking.ranked.tsv")

    # Stage 3 - random region extraction
    stage3 = []
    for seq in stage2:
        n_region = extract_random_region(seq, left_stem, right_stem, n_len)
        if n_region is not None:
            stage3.append(n_region)
    aggregate_and_write(stage3, f"{outdir}/stage3_random_region.ranked.tsv")

    # ------------------------------------------------------------------ #
    # Survival summary
    # ------------------------------------------------------------------ #
    with open(f"{outdir}/survival_summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_reads_processed\t{total_reads}\n")
        out.write(f"raw_unique_sequences\t{raw_unique}\n")
        out.write(f"stage1_matching_survivors\t{len(stage1)}\n")
        out.write(f"stage2_flanking_survivors\t{len(stage2)}\n")
        out.write(f"stage3_random_region_survivors\t{len(stage3)}\n")

    print("AptaSelect complete.", file=sys.stderr)
    print(f"  total reads processed : {total_reads}", file=sys.stderr)
    print(f"  stage1 survivors      : {len(stage1)}", file=sys.stderr)
    print(f"  stage2 survivors      : {len(stage2)}", file=sys.stderr)
    print(f"  stage3 survivors (N)  : {len(stage3)}", file=sys.stderr)


if __name__ == "__main__":
    main()