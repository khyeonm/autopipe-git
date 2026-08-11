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

Optimizations
-------------
OPT 1 (multi-core): each read pair is processed independently, so per-read work
is distributed across worker processes. Results are merged in original read order
so the stable-sort tie ordering (first-seen order) is preserved exactly.

OPT 2 (integer/XOR sequence comparison): each sequence is converted to a single
big integer ONCE (one byte per character); comparing two equal-length regions is
XOR + counting nonzero bytes with bit operations, with sub-regions extracted by
shift/mask. No Python per-character loop, no substring building in the hot loop.

OPT 3 (early-exit overlap search): join searches from the LONGEST overlap down
and stops once the best-attainable score at the current (and every shorter)
length can no longer beat the best found. Tie handling reproduces the original
selection; the per-base merge is skipped when the overlap matches exactly.

OPT 4 (pigeonhole seed-and-verify pattern search): to locate a primer allowing k
mismatches, split it into k+1 non-overlapping pieces. Any window matching within
k mismatches must contain at least one exactly-matching piece, so a C-level exact
search (str.find) over each piece yields a small candidate set; only those
candidates are Hamming-verified (integer/XOR), emitted in the original scan order.

Progress logging is emitted by the PARENT process only, to stderr and to an
optional separate file. It never touches the four TSVs or summary.txt, which stay
byte-for-byte identical to the original.
"""

import argparse
import gzip
import sys
import os
import math
import time
import multiprocessing as mp
from functools import lru_cache
from collections import defaultdict, Counter


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
            name1 = f1.readline().rstrip()
            if not name1:
                break
            seq1 = f1.readline().rstrip()
            f1.readline()
            qual1 = f1.readline().rstrip()
            name2 = f2.readline().rstrip()
            seq2 = f2.readline().rstrip()
            f2.readline()
            qual2 = f2.readline().rstrip()
            yield (name1, seq1, qual1, seq2, qual2)


# ── Progress logging (parent process only; separate from output files) ───────

class _Progress:
    """Timestamped progress lines to stderr (visible live in the run log) and,
    optionally, to a separate file. Writes nothing to the pipeline's output
    files, so results stay byte-for-byte identical."""

    def __init__(self, path="", interval=2.0):
        self.interval = float(interval) if interval else 0.0
        self.start = time.monotonic()
        self.last = self.start
        self.fh = None
        if path:
            try:
                self.fh = open(path, "w", buffering=1)
            except OSError:
                self.fh = None

    def _emit(self, msg):
        line = f"[progress +{time.monotonic() - self.start:8.1f}s] {msg}"
        print(line, file=sys.stderr, flush=True)
        if self.fh:
            self.fh.write(line + "\n")

    def log(self, msg):
        self._emit(msg)
        self.last = time.monotonic()

    def tick(self, msg_fn):
        """Emit at most once per `interval` seconds; `msg_fn` is only evaluated
        when a line is actually written."""
        now = time.monotonic()
        if self.interval <= 0 or (now - self.last) >= self.interval:
            self._emit(msg_fn())
            self.last = now

    def elapsed(self):
        return time.monotonic() - self.start

    def close(self):
        if self.fh:
            try:
                self.fh.close()
            except OSError:
                pass


# ── Integer/XOR sequence comparison (OPT 2) ──────────────────────────────────

def _encode(seq):
    """String -> big-endian integer, one byte per character (latin-1: 1:1 for
    all code points 0..255, covering all FASTQ sequence characters)."""
    return int.from_bytes(seq.encode("latin-1"), "big")


@lru_cache(maxsize=None)
def _field_consts(nbytes):
    """(mask, ones): mask=0xFF..FF, ones=0x0101..01, each nbytes wide."""
    mask = (1 << (8 * nbytes)) - 1
    ones = mask // 0xFF
    return mask, ones


def _count_nonzero_bytes(x, ones):
    """Count nonzero bytes in x. OR-reduce each byte down to its low bit
    (cross-byte contamination only reaches discarded high bits), mask, popcount."""
    y = x
    y |= y >> 4
    y |= y >> 2
    y |= y >> 1
    y &= ones
    return y.bit_count()


@lru_cache(maxsize=None)
def _pattern_consts(pattern):
    pb = pattern.encode("latin-1")
    plen = len(pb)
    pat_int = int.from_bytes(pb, "big")
    mask, ones = _field_consts(plen)
    return pat_int, plen, mask, ones


@lru_cache(maxsize=None)
def _pattern_pieces(pattern, max_mismatches):
    """Split pattern into (max_mismatches + 1) non-overlapping contiguous pieces
    that fully cover it, as evenly as possible. Returns a tuple of (piece,
    offset). Cached because the same primers are used for every read. Only valid
    when len(pattern) >= max_mismatches + 1 (guaranteed by the caller)."""
    plen = len(pattern)
    num_pieces = max_mismatches + 1
    base = plen // num_pieces
    rem = plen % num_pieces
    pieces = []
    offset = 0
    for j in range(num_pieces):
        size = base + (1 if j < rem else 0)
        pieces.append((pattern[offset:offset + size], offset))
        offset += size
    return tuple(pieces)


def hamming_distance(s1, s2):
    """Reference only; hot paths use the integer/XOR method."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


# ── Stage 1: Join ────────────────────────────────────────────────────────────

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

    len_first = len(first_seq)
    len_second = len(second_seq)
    max_overlap = min(len_first, len_second)

    first_int = _encode(first_seq)
    second_int = _encode(second_seq)

    best_score = float('inf')
    best_overlap = -1
    best_d = -1

    for ov in range(max_overlap, min_overlap - 1, -1):
        if 1000.0 / ov > best_score:
            break
        mask, ones = _field_consts(ov)
        tail_int = first_int & mask
        head_int = second_int >> (8 * (len_second - ov))
        d = _count_nonzero_bytes(tail_int ^ head_int, ones)
        mismatch_rate = d / ov
        if mismatch_rate <= max_mismatch_pct:
            score = (1000 * (d * d + 1)) / ov
            if score <= best_score:
                best_score = score
                best_overlap = ov
                best_d = d

    if best_overlap < 0:
        return None

    ov = best_overlap
    first_tail_start = len_first - ov

    if best_d == 0:
        return first_seq + second_seq[ov:]

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


# ── Stages 2–4: Pattern matching ─────────────────────────────────────────────

def pattern_matches(seq, pattern, max_mismatches):
    """Yield start positions where pattern matches seq within `max_mismatches`,
    ascending (identical to a full left-to-right scan). OPT 4 pigeonhole
    seed-and-verify."""
    plen = len(pattern)
    n = len(seq)
    if plen == 0:
        for i in range(n + 1):
            yield i
        return
    if n < plen:
        return
    if max_mismatches >= plen:
        for i in range(n - plen + 1):
            yield i
        return

    last_start = n - plen
    pieces = _pattern_pieces(pattern, max_mismatches)

    candidates = set()
    find = seq.find
    for piece, offset in pieces:
        start = 0
        while True:
            p = find(piece, start)
            if p == -1:
                break
            i = p - offset
            if 0 <= i <= last_start:
                candidates.add(i)
            start = p + 1

    if not candidates:
        return

    seq_int = _encode(seq)
    pat_int, _plen, mask, ones = _pattern_consts(pattern)
    for i in sorted(candidates):
        shift = 8 * (n - i - plen)
        window_int = (seq_int >> shift) & mask
        d = _count_nonzero_bytes(window_int ^ pat_int, ones)
        if d <= max_mismatches:
            yield i


def find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches=1,
                      required_between_length=None):
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
    result = find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches)
    if result is None:
        return None
    left_start, right_end = result
    return seq[left_start:right_end]


def stage3_sort1(seq, left_pattern, right_pattern, max_mismatches):
    result = find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches)
    return seq if result is not None else None


def stage4_sort2(seq, left_pattern, right_pattern, max_mismatches, between_length):
    result = find_pattern_pair(seq, left_pattern, right_pattern, max_mismatches,
                               required_between_length=between_length)
    return seq if result is not None else None


# ── Stage 5: Aggregation & Ranking ───────────────────────────────────────────

def aggregate_and_rank(counter):
    """Sort by count descending. sorted() is stable, so equal-count sequences
    keep insertion (first-seen) order."""
    return sorted(counter.items(), key=lambda x: -x[1])


# ── Core detection (container-aware) ─────────────────────────────────────────

def _cgroup_quota_cores():
    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:
            parts = fh.read().split()
        if parts and parts[0] != "max":
            quota = int(parts[0])
            period = int(parts[1]) if len(parts) > 1 else 100000
            if quota > 0 and period > 0:
                return max(1, int(math.floor(quota / period)))
    except (OSError, ValueError):
        pass
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as fh:
            quota = int(fh.read().strip())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as fh:
            period = int(fh.read().strip())
        if quota > 0 and period > 0:
            return max(1, int(math.floor(quota / period)))
    except (OSError, ValueError):
        pass
    return None


def detect_usable_cores():
    try:
        affinity = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity = os.cpu_count() or 1
    usable = affinity
    quota = _cgroup_quota_cores()
    if quota is not None:
        usable = min(usable, quota)
    return max(1, usable)


def resolve_cores(requested):
    if requested is not None and requested > 0:
        return int(requested)
    return detect_usable_cores()


# ── Worker: process one contiguous block of read pairs ───────────────────────

_PARAMS = None


def _init_worker(params):
    global _PARAMS
    _PARAMS = params


def _process_work_unit(read_tuples):
    p = _PARAMS
    stage1 = {}
    stage2 = {}
    stage3 = {}
    stage4 = {}

    total_pairs = 0
    joined_count = 0
    s2_pass = 0
    s3_pass = 0
    s4_pass = 0

    chunk_size = p["chunk_size"]
    chunk_seqs = []

    sel_left = p["sel_left"]; sel_right = p["sel_right"]
    sort1_left = p["sort1_left"]; sort1_right = p["sort1_right"]
    sort2_left = p["sort2_left"]; sort2_right = p["sort2_right"]
    max_mm = p["max_mismatches"]
    between = p["sort2_between_length"]
    min_overlap = p["min_overlap"]
    max_mismatch_pct = p["max_mismatch_pct"]
    short_mode = p["short_mode"]

    def process_chunk(seqs):
        nonlocal joined_count, s2_pass, s3_pass, s4_pass
        for seq in seqs:
            stage1[seq] = stage1.get(seq, 0) + 1

            extracted = stage2_selection(seq, sel_left, sel_right, max_mm)
            if extracted is None:
                continue
            s2_pass += 1
            stage2[extracted] = stage2.get(extracted, 0) + 1

            passed = stage3_sort1(extracted, sort1_left, sort1_right, max_mm)
            if passed is None:
                continue
            s3_pass += 1
            stage3[passed] = stage3.get(passed, 0) + 1

            passed = stage4_sort2(extracted, sort2_left, sort2_right, max_mm, between)
            if passed is None:
                continue
            s4_pass += 1
            stage4[passed] = stage4.get(passed, 0) + 1

    for (seq1, qual1, seq2, qual2) in read_tuples:
        total_pairs += 1
        joined = join_reads(seq1, qual1, seq2, qual2,
                            min_overlap=min_overlap,
                            max_mismatch_pct=max_mismatch_pct,
                            short_mode=short_mode)
        if joined is None:
            continue
        joined_count += 1
        chunk_seqs.append(joined)
        if len(chunk_seqs) >= chunk_size:
            process_chunk(chunk_seqs)
            chunk_seqs = []
    if chunk_seqs:
        process_chunk(chunk_seqs)

    return ((total_pairs, joined_count, s2_pass, s3_pass, s4_pass),
            stage1, stage2, stage3, stage4)


# ── Main pipeline ────────────────────────────────────────────────────────────

def write_ranked(filepath, ranked):
    with open(filepath, "w") as f:
        f.write("rank\tcount\tsequence\n")
        for rank, (seq, count) in enumerate(ranked, 1):
            f.write(f"{rank}\t{count}\t{seq}\n")


def _iter_work_units(r1_path, r2_path, work_unit_size):
    unit = []
    for rec in read_fastq_pairs(r1_path, r2_path):
        unit.append((rec[1], rec[2], rec[3], rec[4]))
        if len(unit) >= work_unit_size:
            yield unit
            unit = []
    if unit:
        yield unit


def run_pipeline(args):
    prog = _Progress(getattr(args, "progress_log", "") or "",
                     getattr(args, "progress_interval", 2.0))

    n_cores = resolve_cores(args.cores)
    work_unit_size = args.work_unit_size
    if work_unit_size is None or work_unit_size <= 0:
        work_unit_size = 10000

    prog.log(f"run started | r1={os.path.basename(args.r1)} r2={os.path.basename(args.r2)} "
             f"| mode={'short' if args.short_mode else 'long'} | workers={n_cores} "
             f"| work_unit_size={work_unit_size} | chunk_size={args.chunk_size}")

    params = {
        "min_overlap": args.min_overlap,
        "max_mismatch_pct": args.max_mismatch_pct,
        "short_mode": args.short_mode,
        "sel_left": args.sel_left,
        "sel_right": args.sel_right,
        "sort1_left": args.sort1_left,
        "sort1_right": args.sort1_right,
        "sort2_left": args.sort2_left,
        "sort2_right": args.sort2_right,
        "max_mismatches": args.max_mismatches,
        "sort2_between_length": args.sort2_between_length,
        "chunk_size": args.chunk_size,
    }

    stage1_counter = {}
    stage2_counter = {}
    stage3_counter = {}
    stage4_counter = {}

    total_pairs = 0
    joined_count = 0
    s2_pass = 0
    s3_pass = 0
    s4_pass = 0

    def merge(result):
        nonlocal total_pairs, joined_count, s2_pass, s3_pass, s4_pass
        (t, j, p2, p3, p4), c1, c2, c3, c4 = result
        total_pairs += t
        joined_count += j
        s2_pass += p2
        s3_pass += p3
        s4_pass += p4
        for seq, n in c1.items():
            stage1_counter[seq] = stage1_counter.get(seq, 0) + n
        for seq, n in c2.items():
            stage2_counter[seq] = stage2_counter.get(seq, 0) + n
        for seq, n in c3.items():
            stage3_counter[seq] = stage3_counter.get(seq, 0) + n
        for seq, n in c4.items():
            stage4_counter[seq] = stage4_counter.get(seq, 0) + n

    def progress_msg(units_done):
        el = prog.elapsed()
        rate = total_pairs / el if el > 0 else 0.0
        return (f"processed {units_done} work units | {total_pairs:,} read pairs "
                f"| joined {joined_count:,} | selection {s2_pass:,} | sort1 {s3_pass:,} "
                f"| sort2 {s4_pass:,} | {rate:,.0f} pairs/s")

    print(f"Processing {args.r1} and {args.r2} ...")
    print(f"Short mode: {args.short_mode}")
    print(f"Worker processes: {n_cores}  (work_unit_size={work_unit_size} read pairs)")

    units = _iter_work_units(args.r1, args.r2, work_unit_size)
    units_done = 0

    if n_cores <= 1:
        _init_worker(params)
        for unit in units:
            merge(_process_work_unit(unit))
            units_done += 1
            if units_done == 1:
                prog.log(progress_msg(units_done))
            else:
                prog.tick(lambda: progress_msg(units_done))
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=n_cores,
                      initializer=_init_worker,
                      initargs=(params,)) as pool:
            for result in pool.imap(_process_work_unit, units, chunksize=1):
                merge(result)
                units_done += 1
                if units_done == 1:
                    prog.log(progress_msg(units_done))
                else:
                    prog.tick(lambda: progress_msg(units_done))

    prog.log(f"read phase complete | {units_done} work units | {total_pairs:,} read pairs "
             f"| joined {joined_count:,} | aggregating and writing outputs ...")

    print(f"\n=== AptaSelect Summary ===")
    print(f"Total read pairs:           {total_pairs}")
    print(f"Stage 1 (Join):             {joined_count} ({joined_count/max(total_pairs,1)*100:.1f}%)")
    print(f"Stage 2 (Selection Filter): {s2_pass} ({s2_pass/max(total_pairs,1)*100:.1f}%)")
    print(f"Stage 3 (1st Sort Filter):  {s3_pass} ({s3_pass/max(total_pairs,1)*100:.1f}%)")
    print(f"Stage 4 (2nd Sort Filter):  {s4_pass} ({s4_pass/max(total_pairs,1)*100:.1f}%)")

    os.makedirs(args.outdir, exist_ok=True)

    ranked1 = aggregate_and_rank(stage1_counter)
    prog.log(f"writing stage1_joined_ranked.tsv ({len(ranked1):,} unique sequences)")
    write_ranked(os.path.join(args.outdir, "stage1_joined_ranked.tsv"), ranked1)
    print(f"Stage 1 unique sequences:   {len(ranked1)}")

    ranked2 = aggregate_and_rank(stage2_counter)
    prog.log(f"writing stage2_selection_ranked.tsv ({len(ranked2):,} unique sequences)")
    write_ranked(os.path.join(args.outdir, "stage2_selection_ranked.tsv"), ranked2)
    print(f"Stage 2 unique sequences:   {len(ranked2)}")

    ranked3 = aggregate_and_rank(stage3_counter)
    prog.log(f"writing stage3_sort1_ranked.tsv ({len(ranked3):,} unique sequences)")
    write_ranked(os.path.join(args.outdir, "stage3_sort1_ranked.tsv"), ranked3)
    print(f"Stage 3 unique sequences:   {len(ranked3)}")

    ranked4 = aggregate_and_rank(stage4_counter)
    prog.log(f"writing stage4_sort2_ranked.tsv ({len(ranked4):,} unique sequences)")
    write_ranked(os.path.join(args.outdir, "stage4_sort2_ranked.tsv"), ranked4)
    print(f"Stage 4 unique sequences:   {len(ranked4)}")

    prog.log("writing summary.txt")
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

    prog.log(f"done | {total_pairs:,} read pairs | joined {joined_count:,} "
             f"| unique s1={len(ranked1):,} s2={len(ranked2):,} s3={len(ranked3):,} "
             f"s4={len(ranked4):,} | total {prog.elapsed():.1f}s")
    prog.close()

    print("\nDone! Results written to:", args.outdir)


def main():
    parser = argparse.ArgumentParser(description="AptaSelect: SELEX aptamer candidate identification")

    parser.add_argument("--r1", required=True, help="Read 1 FASTQ file path")
    parser.add_argument("--r2", required=True, help="Read 2 FASTQ file path")
    parser.add_argument("--outdir", required=True, help="Output directory")

    parser.add_argument("--short-mode", action="store_true", default=True,
                        help="Use short library mode (default: True)")
    parser.add_argument("--long-mode", action="store_true", default=False,
                        help="Use long library mode (overrides short-mode)")
    parser.add_argument("--min-overlap", type=int, default=6,
                        help="Minimum overlap length in bp (default: 6)")
    parser.add_argument("--max-mismatch-pct", type=float, default=0.08,
                        help="Maximum mismatch rate for overlap (default: 0.08)")

    parser.add_argument("--max-mismatches", type=int, default=1,
                        help="Maximum mismatches for pattern matching (default: 1)")

    parser.add_argument("--sel-left", default="CCACTTCTCCTTCCATCCTAAAC",
                        help="Selection left pattern")
    parser.add_argument("--sel-right", default="GAGTAGTTTGGAGGGTTGTCTG",
                        help="Selection right pattern")

    parser.add_argument("--sort1-left", default="TCCTAAAC",
                        help="1st Sort left pattern")
    parser.add_argument("--sort1-right", default="GAGTAGTT",
                        help="1st Sort right pattern")

    parser.add_argument("--sort2-left", default="TCTCTCTCTC",
                        help="2nd Sort left pattern")
    parser.add_argument("--sort2-right", default="GAGAGAGAGA",
                        help="2nd Sort right pattern")
    parser.add_argument("--sort2-between-length", type=int, default=20,
                        help="Required between-length for 2nd Sort (default: 20)")

    parser.add_argument("--chunk-size", type=int, default=10000,
                        help="Processing chunk size (default: 10000)")

    parser.add_argument("--cores", type=int, default=0,
                        help="Worker processes. 0 or negative = auto-detect (default: 0)")
    parser.add_argument("--work-unit-size", type=int, default=10000,
                        help="Read pairs per worker task; separate from --chunk-size (default: 10000)")

    parser.add_argument("--progress-log", default="",
                        help="Optional path for a separate progress log file. "
                             "Progress is always also written to stderr.")
    parser.add_argument("--progress-interval", type=float, default=2.0,
                        help="Minimum seconds between progress lines (default: 2.0)")

    args = parser.parse_args()

    if args.long_mode:
        args.short_mode = False

    run_pipeline(args)


if __name__ == "__main__":
    main()