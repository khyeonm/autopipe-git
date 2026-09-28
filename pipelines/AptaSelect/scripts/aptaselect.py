#!/usr/bin/env python3
"""AptaSelect baseline: single-process, pure-Python aptamer candidate selection.

Each input read (Read 1 and Read 2 alike) is handled as its own sequence and
passed through three sequential stages:

  Stage 1  matching         - LCR and RCR both present (edit distance <= k),
                              read or its reverse complement.
  Stage 2  flanking         - exact left/right markers, extract what lies between.
  Stage 3  random region    - exact left/right stems with exactly N bases between.

Every stage (plus the raw input as "stage 0") is aggregated into a ranked table
(sequence<TAB>count, most frequent first, ties kept in first-seen order).

This is deliberately the unoptimised baseline: no multiprocessing and no
compiled / external string-matching libraries.
"""

import argparse
import gzip
import os
import sys
import time

ORDINARY = frozenset("ACGT")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# --------------------------------------------------------------------------- IO
def open_text(path):
    """Open a FASTQ file, transparently handling gzip (detected by magic bytes)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fastq_seqs(path):
    """Yield the sequence line of every FASTQ record (upper-cased)."""
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            if not header.strip():
                continue
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not header.startswith("@") or not plus.startswith("+") or not qual:
                sys.exit(f"ERROR: malformed FASTQ record in {path}: {header.strip()!r}")
            yield seq.strip().upper()


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------- Stage 1 logic
def approx_contains(pattern, text, max_err):
    """True if the whole `pattern` occurs somewhere in `text` with edit distance
    (substitution, insertion or deletion, each costing 1) <= max_err.

    Plain-Python semi-global dynamic programming (Sellers' algorithm): the
    pattern must be aligned end to end, the matched stretch of text may start
    and end anywhere. Column-wise over the text, O(len(pattern) * len(text)).
    """
    m = len(pattern)
    if m == 0:
        return True
    if max_err >= m:
        return True  # the pattern can always be deleted entirely
    if max_err == 0:
        return pattern in text
    prev = list(range(m + 1))  # column before any text character
    for ch in text:
        cur = [0] * (m + 1)  # row 0 stays 0: alignment may start at any text position
        for i in range(1, m + 1):
            cost = 0 if pattern[i - 1] == ch else 1
            best = prev[i - 1] + cost          # match / substitution
            ins = prev[i] + 1                  # extra base in text
            if ins < best:
                best = ins
            dele = cur[i - 1] + 1              # base of pattern missing in text
            if dele < best:
                best = dele
            cur[i] = best
        if cur[m] <= max_err:
            return True
        prev = cur
    return False


def stage1_match(read, lcr, rcr, max_err):
    """Return the read (or its reverse complement) oriented LCR...RCR, else None."""
    if approx_contains(lcr, read, max_err) and approx_contains(rcr, read, max_err):
        return read
    rc = revcomp(read)
    if approx_contains(lcr, rc, max_err) and approx_contains(rcr, rc, max_err):
        return rc
    return None


# ---------------------------------------------------------------- Stage 2 logic
def find_all(text, pattern):
    """All (overlapping) start positions of an exact `pattern` in `text`."""
    positions = []
    start = text.find(pattern)
    while start != -1:
        positions.append(start)
        start = text.find(pattern, start + 1)
    return positions


def stage2_flank(seq, left_marker, right_marker):
    """Extract the sequence between the left and right markers (exact match).

    Try left-marker occurrences from leftmost onward. From the end of a left
    marker, the uninterrupted run of ordinary bases (A/C/G/T) extends up to the
    first other character (e.g. N). Among right markers lying entirely inside
    that run, take the farthest one. If none is reachable, move on to the next
    left-marker occurrence. No length limit, no stem check.
    """
    if left_marker not in seq or right_marker not in seq:
        return None
    right_positions = find_all(seq, right_marker)
    n = len(seq)
    rlen = len(right_marker)
    for lpos in find_all(seq, left_marker):
        start = lpos + len(left_marker)
        run_end = start
        while run_end < n and seq[run_end] in ORDINARY:
            run_end += 1
        best = None
        for rpos in right_positions:
            if rpos >= start and rpos + rlen <= run_end:
                best = rpos  # positions are ascending -> keeps the farthest
        if best is not None:
            return seq[start:best]
    return None


# ---------------------------------------------------------------- Stage 3 logic
def stage3_random_region(seq, left_stem, right_stem, n_len):
    """Return the exactly-n_len region between the stems, else None.

    All left/right stem position pairs are considered (overlapping matches
    included); pairs are visited left-stem position first, then right-stem
    position, both ascending, and the first pair whose gap equals n_len wins.
    """
    llen = len(left_stem)
    rlen = len(right_stem)
    right_set = set(find_all(seq, right_stem))
    if not right_set:
        return None
    for lpos in find_all(seq, left_stem):
        n_start = lpos + llen
        rpos = n_start + n_len
        if rpos + rlen <= len(seq) and rpos in right_set:
            return seq[n_start:rpos]
    return None


# ----------------------------------------------------------------- aggregation
class OrderedCounter:
    """Counts occurrences while remembering first-seen order (dicts keep order)."""

    def __init__(self):
        self.counts = {}
        self.total = 0

    def add(self, seq):
        self.counts[seq] = self.counts.get(seq, 0) + 1
        self.total += 1

    def ranked(self):
        # sorted() is stable: equal counts stay in first-seen order.
        return sorted(self.counts.items(), key=lambda kv: -kv[1])


def write_ranked(counter, path):
    with open(path, "w") as out:
        out.write("sequence\tcount\n")
        for seq, count in counter.ranked():
            out.write(f"{seq}\t{count}\n")


# ------------------------------------------------------------------------ main
def check_seq(name, value, allow_empty=False):
    value = (value or "").strip().upper()
    if not value and not allow_empty:
        sys.exit(f"ERROR: config value '{name}' is empty; please set it.")
    bad = set(value) - ORDINARY
    if bad:
        sys.exit(f"ERROR: config value '{name}' contains non-ACGT characters: {sorted(bad)}")
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reads", nargs="+", required=True, help="FASTQ(.gz) files; each read is its own sequence")
    p.add_argument("--outdir", required=True)
    p.add_argument("--lcr", required=True)
    p.add_argument("--rcr", required=True)
    p.add_argument("--max-errors", type=int, required=True)
    p.add_argument("--left-marker", required=True)
    p.add_argument("--right-marker", required=True)
    p.add_argument("--left-stem", required=True)
    p.add_argument("--right-stem", required=True)
    p.add_argument("--random-length", type=int, required=True)
    p.add_argument("--progress-every", type=int, default=100000)
    a = p.parse_args()

    lcr = check_seq("lcr", a.lcr)
    rcr = check_seq("rcr", a.rcr)
    left_marker = check_seq("left_marker", a.left_marker)
    right_marker = check_seq("right_marker", a.right_marker)
    left_stem = check_seq("left_stem", a.left_stem)
    right_stem = check_seq("right_stem", a.right_stem)
    if a.max_errors < 0:
        sys.exit("ERROR: max_errors must be >= 0")
    if a.random_length < 1:
        sys.exit("ERROR: random_region_length must be >= 1")

    insert_len = len(lcr) + len(left_stem) + a.random_length + len(right_stem) + len(rcr)

    os.makedirs(a.outdir, exist_ok=True)
    raw, s1, s2, s3 = OrderedCounter(), OrderedCounter(), OrderedCounter(), OrderedCounter()

    t0 = time.time()
    processed = 0
    for path in a.reads:
        print(f"[aptaselect] reading {path}", file=sys.stderr, flush=True)
        for read in iter_fastq_seqs(path):
            processed += 1
            raw.add(read)
            oriented = stage1_match(read, lcr, rcr, a.max_errors)
            if oriented is not None:
                s1.add(oriented)
                flank = stage2_flank(oriented, left_marker, right_marker)
                if flank is not None:
                    s2.add(flank)
                    n_region = stage3_random_region(flank, left_stem, right_stem, a.random_length)
                    if n_region is not None:
                        s3.add(n_region)
            if a.progress_every and processed % a.progress_every == 0:
                print(f"[aptaselect] {processed} reads, {time.time() - t0:.1f}s", file=sys.stderr, flush=True)

    write_ranked(raw, os.path.join(a.outdir, "stage0_raw.ranked.tsv"))
    write_ranked(s1, os.path.join(a.outdir, "stage1_matching.ranked.tsv"))
    write_ranked(s2, os.path.join(a.outdir, "stage2_flanking.ranked.tsv"))
    write_ranked(s3, os.path.join(a.outdir, "stage3_random_region.ranked.tsv"))

    summary = [
        ("total_reads_processed", processed),
        ("raw_unique_sequences", len(raw.counts)),
        ("stage1_matching_survivors", s1.total),
        ("stage2_flanking_survivors", s2.total),
        ("stage3_random_region_survivors", s3.total),
        ("derived_insert_length", insert_len),
    ]
    with open(os.path.join(a.outdir, "survival_summary.tsv"), "w") as out:
        for name, value in summary:
            out.write(f"{name}\t{value}\n")

    for name, value in summary:
        print(f"[aptaselect] {name}\t{value}", file=sys.stderr)
    print(f"[aptaselect] done in {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()