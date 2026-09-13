#!/usr/bin/env python3
"""AptaSelect: identify high-frequency aptamer candidate sequences from paired-end
SELEX FASTQ files via three sequential filtering stages plus aggregation.

Every pattern sequence, mismatch tolerance and length is passed in from config.yaml
(via command-line arguments) — nothing about the library design is hard-coded here.

Insert layout: LCR - left_stem - N(random) - right_stem - RCR
Stem lengths follow from their sequences; the overall insert length is derived,
never hard-coded.

R1 and R2 are each treated as their own sequence, processed independently through
all three stages, and their survivors are pooled into a single per-stage result set.
"""

import argparse
import gzip
import sys
from collections import OrderedDict


# ----------------------------------------------------------------------------- #
# Basic sequence utilities
# ----------------------------------------------------------------------------- #

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    """Return the reverse complement of a nucleotide sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def open_maybe_gzip(path):
    """Open a file, transparently handling gzip by magic number."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fastq(path):
    """Yield sequence strings (uppercased) from a FASTQ file."""
    with open_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                break  # truncated record
            yield seq.strip().upper()


# ----------------------------------------------------------------------------- #
# Approximate matching (Hamming distance within a mismatch tolerance)
# ----------------------------------------------------------------------------- #

def hamming_within(text, pattern, start, max_mismatch):
    """True if pattern matches text starting at `start` within max_mismatch
    substitutions (equal-length window; no indels)."""
    end = start + len(pattern)
    if end > len(text):
        return False
    mism = 0
    for a, b in zip(text[start:end], pattern):
        if a != b:
            mism += 1
            if mism > max_mismatch:
                return False
    return True


def find_approx(text, pattern, max_mismatch, from_pos=0):
    """Return the leftmost start index >= from_pos where pattern matches text
    within max_mismatch substitutions, or -1 if none."""
    plen = len(pattern)
    last_start = len(text) - plen
    for start in range(from_pos, last_start + 1):
        if hamming_within(text, pattern, start, max_mismatch):
            return start
    return -1


def find_exact_all(text, pattern, from_pos=0):
    """Yield every start index (>= from_pos) of an exact (possibly overlapping)
    occurrence of pattern in text, left to right."""
    plen = len(pattern)
    if plen == 0:
        return
    start = from_pos
    while True:
        idx = text.find(pattern, start)
        if idx == -1:
            return
        yield idx
        start = idx + 1  # allow overlapping matches


# ----------------------------------------------------------------------------- #
# Stage 1 - Matching: orient by LCR ... RCR within mismatch tolerance
# ----------------------------------------------------------------------------- #

def stage1_orient(seq, lcr, rcr, mismatch):
    """Return the oriented sequence if it contains LCR followed downstream by RCR
    (each within `mismatch` substitutions), checking the read as-is first and then
    its reverse complement. Return None if neither orientation qualifies."""
    for candidate in (seq, reverse_complement(seq)):
        lcr_pos = find_approx(candidate, lcr, mismatch, from_pos=0)
        if lcr_pos == -1:
            continue
        # RCR must lie downstream of (after) the LCR match
        rcr_search_from = lcr_pos + len(lcr)
        rcr_pos = find_approx(candidate, rcr, mismatch, from_pos=rcr_search_from)
        if rcr_pos != -1:
            return candidate
    return None


# ----------------------------------------------------------------------------- #
# Stage 2 - Flanking extraction: between two exact-match markers
# ----------------------------------------------------------------------------- #

def stage2_extract(seq, left_marker, right_marker):
    """Find left_marker (inner end of LCR) and right_marker (inner start of RCR)
    by exact match, and return the sequence lying strictly between them. The
    right marker is searched downstream of the end of the left marker. Length is
    NOT checked here. Return None if either marker is not found exactly."""
    lpos = seq.find(left_marker)
    if lpos == -1:
        return None
    inner_start = lpos + len(left_marker)
    rpos = seq.find(right_marker, inner_start)
    if rpos == -1:
        return None
    return seq[inner_start:rpos]


# ----------------------------------------------------------------------------- #
# Stage 3 - N (random region) extraction: exact stems, exact between-length
# ----------------------------------------------------------------------------- #

def stage3_random_region(extracted, left_stem, right_stem, target_len):
    """Find left_stem and right_stem by exact match and require the region between
    them to be exactly target_len. Consider all overlapping left/right position
    pairs and accept the FIRST pair (left ascending, then right ascending) that
    satisfies the exact-length constraint. Return the N region, or None."""
    left_positions = list(find_exact_all(extracted, left_stem, from_pos=0))
    if not left_positions:
        return None
    llen = len(left_stem)
    for lpos in left_positions:
        n_start = lpos + llen
        # right stem must begin at or after the end of the left stem
        for rpos in find_exact_all(extracted, right_stem, from_pos=n_start):
            if rpos - n_start == target_len:
                return extracted[n_start:rpos]
    return None


# ----------------------------------------------------------------------------- #
# Aggregation and ranking
# ----------------------------------------------------------------------------- #

def aggregate_ranked(sequences):
    """Count identical sequences and return a list of (sequence, count) sorted by
    descending count, ties broken by first-seen order (stable). No secondary
    alphabetical key is applied."""
    counts = OrderedDict()  # preserves first-seen order
    for s in sequences:
        counts[s] = counts.get(s, 0) + 1
    items = list(counts.items())  # already in first-seen order
    items.sort(key=lambda kv: kv[1], reverse=True)  # stable sort on count desc
    return items


def write_tsv(path, ranked_items):
    with open(path, "w") as out:
        out.write("sequence\tcount\n")
        for seq, count in ranked_items:
            out.write("{}\t{}\n".format(seq, count))


# ----------------------------------------------------------------------------- #
# Driver
# ----------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="AptaSelect three-stage SELEX filter")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--lcr", required=True)
    ap.add_argument("--rcr", required=True)
    ap.add_argument("--lcr-rcr-mismatch", type=int, required=True)
    ap.add_argument("--left-marker", required=True)
    ap.add_argument("--right-marker", required=True)
    ap.add_argument("--left-stem", required=True)
    ap.add_argument("--right-stem", required=True)
    ap.add_argument("--random-region-length", type=int, required=True)
    ap.add_argument("--out-stage1", required=True)
    ap.add_argument("--out-stage2", required=True)
    ap.add_argument("--out-stage3", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    lcr = args.lcr.upper()
    rcr = args.rcr.upper()
    left_marker = args.left_marker.upper()
    right_marker = args.right_marker.upper()
    left_stem = args.left_stem.upper()
    right_stem = args.right_stem.upper()
    mismatch = args.lcr_rcr_mismatch
    target_len = args.random_region_length

    # Derived (informational only): full insert length is never hard-coded.
    derived_insert_len = (
        len(lcr) + len(left_stem) + target_len + len(right_stem) + len(rcr)
    )

    stage1_seqs = []  # oriented reads passing Stage 1 (R1 and R2 pooled)
    stage2_seqs = []  # extracted regions passing Stage 2
    stage3_seqs = []  # N regions passing Stage 3

    total_reads = 0
    # Both mates handled independently, then pooled into one set per stage.
    for path in (args.r1, args.r2):
        for seq in read_fastq(path):
            total_reads += 1

            oriented = stage1_orient(seq, lcr, rcr, mismatch)
            if oriented is None:
                continue
            stage1_seqs.append(oriented)

            extracted = stage2_extract(oriented, left_marker, right_marker)
            if extracted is None:
                continue
            stage2_seqs.append(extracted)

            n_region = stage3_random_region(
                extracted, left_stem, right_stem, target_len
            )
            if n_region is None:
                continue
            stage3_seqs.append(n_region)

    write_tsv(args.out_stage1, aggregate_ranked(stage1_seqs))
    write_tsv(args.out_stage2, aggregate_ranked(stage2_seqs))
    write_tsv(args.out_stage3, aggregate_ranked(stage3_seqs))

    with open(args.summary, "w") as s:
        s.write("metric\tvalue\n")
        s.write("total_reads_processed\t{}\n".format(total_reads))
        s.write("stage1_matching_survivors\t{}\n".format(len(stage1_seqs)))
        s.write("stage2_flanking_survivors\t{}\n".format(len(stage2_seqs)))
        s.write("stage3_random_region_survivors\t{}\n".format(len(stage3_seqs)))
        s.write("derived_insert_length\t{}\n".format(derived_insert_len))

    # Also echo the survival counts to stderr for the Snakemake log.
    sys.stderr.write(
        "AptaSelect done. reads={} stage1={} stage2={} stage3={} "
        "(derived insert length={})\n".format(
            total_reads,
            len(stage1_seqs),
            len(stage2_seqs),
            len(stage3_seqs),
            derived_insert_len,
        )
    )


if __name__ == "__main__":
    main()