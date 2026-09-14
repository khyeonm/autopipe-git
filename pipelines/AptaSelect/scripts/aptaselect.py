#!/usr/bin/env python3
"""
AptaSelect: identify high-frequency aptamer candidate sequences from paired-end
FASTQ files produced by SELEX experiments.

Stages (each read mate processed independently, survivors pooled):
  Stage 0 - Raw: aggregate all raw reads by identical sequence.
  Stage 1 - Matching: orient each read by LCR...RCR (mismatch-tolerant), both strands.
  Stage 2 - Flanking extraction: region between inner LCR/RCR markers (exact).
  Stage 3 - N extraction: region between stems, exact target length.

Parallelism: reads are split into contiguous, index-ordered tasks and processed
across worker processes; results are reassembled strictly in task order, so
sequences enter every counter in the same order as a single-process run,
preserving the first-seen tie-break of the stable ranking sort -> byte-identical
outputs regardless of worker count or task size.

Fuzzy matching (seed-and-verify / pigeonhole): to find a pattern within k
mismatches, split it into k+1 non-overlapping pieces. Any alignment with <= k
mismatches must match at least one piece exactly (pigeonhole). Each piece is
located with str.find (C-implemented exact search); each hit implies a candidate
start = hit_pos - piece_offset. Candidate starts are gathered, sorted ascending
and de-duplicated, then the true mismatch count is verified only at those
candidates. Returning the first verified candidate reproduces exactly the
leftmost match the naive left-to-right scan would have found. The k+1 split is
cached per (pattern, k) since the same patterns are reused across all reads.

Progress logging: human-readable, timestamped progress is written to STDERR
only. It never touches the five output TSVs, whose bytes are independent of
logging. In the Snakemake rule stderr is redirected to a separate log file.
"""

import argparse
import gzip
import os
import sys
import time
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor


_T0 = time.monotonic()


def log(msg):
    """Timestamped progress line to stderr (separate from output files)."""
    elapsed = time.monotonic() - _T0
    print(f"[aptaselect +{elapsed:7.1f}s] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# FASTQ reading
# --------------------------------------------------------------------------- #
def open_maybe_gzip(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fastq_sequences(path):
    with open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
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


def find_all_exact(text, pattern):
    positions = []
    if not pattern:
        return positions
    start = text.find(pattern)
    while start != -1:
        positions.append(start)
        start = text.find(pattern, start + 1)
    return positions


# --------------------------------------------------------------------------- #
# Seed-and-verify fuzzy search (pigeonhole)
# --------------------------------------------------------------------------- #
_SPLIT_CACHE = {}


def _split_pieces(pattern, k):
    """
    Split pattern into k+1 non-overlapping, contiguous, covering pieces.
    Returns list of (offset, piece). Sizes are as even as possible. Cached.
    """
    key = (pattern, k)
    cached = _SPLIT_CACHE.get(key)
    if cached is not None:
        return cached
    plen = len(pattern)
    n = k + 1
    if plen == 0:
        pieces = [(0, "")]
        _SPLIT_CACHE[key] = pieces
        return pieces
    if n > plen:
        n = plen
    base = plen // n
    rem = plen % n
    pieces = []
    off = 0
    for i in range(n):
        size = base + (1 if i < rem else 0)
        pieces.append((off, pattern[off:off + size]))
        off += size
    _SPLIT_CACHE[key] = pieces
    return pieces


def find_first_fuzzy(text, pattern, max_mismatch, from_pos=0):
    """
    Leftmost start >= from_pos at which pattern matches text within
    max_mismatch substitutions, else -1. Uses seed-and-verify; result is
    identical to a naive left-to-right scan.
    """
    plen = len(pattern)
    tlen = len(text)
    if plen == 0:
        return from_pos
    last_start = tlen - plen
    if from_pos > last_start:
        return -1
    if from_pos < 0:
        from_pos = 0

    if max_mismatch == 0:
        pos = text.find(pattern, from_pos)
        return pos if pos != -1 and pos <= last_start else -1

    pieces = _split_pieces(pattern, max_mismatch)

    candidates = set()
    for off, piece in pieces:
        search_from = from_pos + off
        if search_from < 0:
            search_from = 0
        pos = text.find(piece, search_from)
        while pos != -1:
            start = pos - off
            if from_pos <= start <= last_start:
                candidates.add(start)
            pos = text.find(piece, pos + 1)

    if not candidates:
        return -1

    for start in sorted(candidates):
        if hamming_at_most(text, pattern, start, max_mismatch):
            return start
    return -1


# --------------------------------------------------------------------------- #
# Stage functions
# --------------------------------------------------------------------------- #
def orient_read(seq, lcr, rcr, max_mismatch):
    for candidate in (seq, reverse_complement(seq)):
        lcr_start = find_first_fuzzy(candidate, lcr, max_mismatch, 0)
        if lcr_start == -1:
            continue
        lcr_end = lcr_start + len(lcr)
        rcr_start = find_first_fuzzy(candidate, rcr, max_mismatch, lcr_end)
        if rcr_start != -1:
            return candidate
    return None


def extract_flanking(seq, left_marker, right_marker):
    lpos = seq.find(left_marker)
    if lpos == -1:
        return None
    inner_start = lpos + len(left_marker)
    rpos = seq.find(right_marker, inner_start)
    if rpos == -1:
        return None
    return seq[inner_start:rpos]


def extract_random_region(seq, left_stem, right_stem, target_len):
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
# Worker
# --------------------------------------------------------------------------- #
_PARAMS = {}


def _init_worker(lcr, rcr, const_mm, left_marker, right_marker,
                 left_stem, right_stem, n_len):
    _PARAMS.update(
        lcr=lcr, rcr=rcr, const_mm=const_mm,
        left_marker=left_marker, right_marker=right_marker,
        left_stem=left_stem, right_stem=right_stem, n_len=n_len,
    )
    _split_pieces(lcr, const_mm)
    _split_pieces(rcr, const_mm)


def _process_task(reads):
    lcr = _PARAMS["lcr"]; rcr = _PARAMS["rcr"]; const_mm = _PARAMS["const_mm"]
    left_marker = _PARAMS["left_marker"]; right_marker = _PARAMS["right_marker"]
    left_stem = _PARAMS["left_stem"]; right_stem = _PARAMS["right_stem"]
    n_len = _PARAMS["n_len"]

    stage1 = []; stage2 = []; stage3 = []
    for seq in reads:
        oriented = orient_read(seq, lcr, rcr, const_mm)
        if oriented is None:
            continue
        stage1.append(oriented)
        region = extract_flanking(oriented, left_marker, right_marker)
        if region is None:
            continue
        stage2.append(region)
        n_region = extract_random_region(region, left_stem, right_stem, n_len)
        if n_region is not None:
            stage3.append(n_region)
    return stage1, stage2, stage3


def _iter_tasks(reads, task_size):
    for start in range(0, len(reads), task_size):
        yield reads[start:start + task_size]


# --------------------------------------------------------------------------- #
def detect_usable_cores():
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def aggregate_and_write(sequences, out_path):
    counts = OrderedDict()
    total = 0
    for s in sequences:
        total += 1
        counts[s] = counts.get(s, 0) + 1
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    with open(out_path, "w") as out:
        out.write("sequence\tcount\n")
        for seq, cnt in items:
            out.write(f"{seq}\t{cnt}\n")
    return total, len(counts)


def main():
    ap = argparse.ArgumentParser(description="AptaSelect SELEX aptamer filter")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--lcr", required=True)
    ap.add_argument("--rcr", required=True)
    ap.add_argument("--constant-mismatch", type=int, required=True)
    ap.add_argument("--left-marker", required=True)
    ap.add_argument("--right-marker", required=True)
    ap.add_argument("--left-stem", required=True)
    ap.add_argument("--right-stem", required=True)
    ap.add_argument("--random-region-length", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--task-size", type=int, default=10000)
    args = ap.parse_args()

    outdir = args.outdir.rstrip("/")
    lcr = args.lcr.strip().upper()
    rcr = args.rcr.strip().upper()
    left_marker = args.left_marker.strip().upper()
    right_marker = args.right_marker.strip().upper()
    left_stem = args.left_stem.strip().upper()
    right_stem = args.right_stem.strip().upper()
    const_mm = args.constant_mismatch
    n_len = args.random_region_length

    n_jobs = args.jobs if args.jobs and args.jobs > 0 else detect_usable_cores()
    task_size = max(1, args.task_size)

    log(f"start: jobs={n_jobs}, task_size={task_size}, const_mismatch={const_mm}")

    # ---- Read input ---- #
    log("reading input FASTQ files ...")
    raw_reads = []
    for path in (args.r1, args.r2):
        before = len(raw_reads)
        for seq in read_fastq_sequences(path):
            raw_reads.append(seq.upper())
        log(f"  loaded {len(raw_reads) - before:,} reads from {os.path.basename(path)}")
    n_reads = len(raw_reads)
    log(f"total reads loaded: {n_reads:,}")

    # ---- Stage 0 ---- #
    log("stage 0: aggregating raw reads ...")
    total_reads, raw_unique = aggregate_and_write(
        raw_reads, f"{outdir}/stage0_raw.ranked.tsv"
    )
    log(f"stage 0 done: {total_reads:,} reads, {raw_unique:,} unique")

    # ---- Stages 1-3 (parallel over ordered tasks) ---- #
    n_tasks = (n_reads + task_size - 1) // task_size if n_reads else 0
    log(f"stages 1-3: processing {n_reads:,} reads in {n_tasks:,} task(s) "
        f"across {n_jobs} worker(s) ...")

    stage1 = []; stage2 = []; stage3 = []
    done = 0
    report_every = max(1, n_tasks // 100)  # ~1% granularity

    run_parallel = n_jobs > 1 and n_reads > task_size
    if run_parallel:
        with ProcessPoolExecutor(
            max_workers=n_jobs,
            initializer=_init_worker,
            initargs=(lcr, rcr, const_mm, left_marker, right_marker,
                      left_stem, right_stem, n_len),
        ) as ex:
            for s1, s2, s3 in ex.map(_process_task,
                                     _iter_tasks(raw_reads, task_size)):
                stage1.extend(s1); stage2.extend(s2); stage3.extend(s3)
                done += 1
                if done % report_every == 0 or done == n_tasks:
                    pct = 100.0 * done / n_tasks
                    log(f"  progress: {done:,}/{n_tasks:,} tasks "
                        f"({pct:5.1f}%) | stage1={len(stage1):,} "
                        f"stage2={len(stage2):,} stage3={len(stage3):,}")
    else:
        _init_worker(lcr, rcr, const_mm, left_marker, right_marker,
                     left_stem, right_stem, n_len)
        for task in _iter_tasks(raw_reads, task_size):
            s1, s2, s3 = _process_task(task)
            stage1.extend(s1); stage2.extend(s2); stage3.extend(s3)
            done += 1
            if done % report_every == 0 or done == n_tasks:
                pct = 100.0 * done / n_tasks if n_tasks else 100.0
                log(f"  progress: {done:,}/{n_tasks:,} tasks "
                    f"({pct:5.1f}%) | stage1={len(stage1):,} "
                    f"stage2={len(stage2):,} stage3={len(stage3):,}")

    log(f"stages 1-3 matching done: stage1={len(stage1):,} "
        f"stage2={len(stage2):,} stage3={len(stage3):,}")

    # ---- Write ranked TSVs ---- #
    log("writing stage 1 ranked TSV ...")
    aggregate_and_write(stage1, f"{outdir}/stage1_matching.ranked.tsv")
    log("writing stage 2 ranked TSV ...")
    aggregate_and_write(stage2, f"{outdir}/stage2_flanking.ranked.tsv")
    log("writing stage 3 ranked TSV ...")
    aggregate_and_write(stage3, f"{outdir}/stage3_random_region.ranked.tsv")

    with open(f"{outdir}/survival_summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_reads_processed\t{total_reads}\n")
        out.write(f"raw_unique_sequences\t{raw_unique}\n")
        out.write(f"stage1_matching_survivors\t{len(stage1)}\n")
        out.write(f"stage2_flanking_survivors\t{len(stage2)}\n")
        out.write(f"stage3_random_region_survivors\t{len(stage3)}\n")

    log("complete.")
    print("AptaSelect complete.", file=sys.stderr)
    print(f"  workers used          : {n_jobs}", file=sys.stderr)
    print(f"  task size             : {task_size}", file=sys.stderr)
    print(f"  total reads processed : {total_reads}", file=sys.stderr)
    print(f"  stage1 survivors      : {len(stage1)}", file=sys.stderr)
    print(f"  stage2 survivors      : {len(stage2)}", file=sys.stderr)
    print(f"  stage3 survivors (N)  : {len(stage3)}", file=sys.stderr)


if __name__ == "__main__":
    main()