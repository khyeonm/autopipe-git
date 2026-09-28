#!/usr/bin/env python3
"""AptaSelect: aptamer candidate selection, multi-process version.

Each input read (Read 1 and Read 2 alike) is handled as its own sequence and
passed through three sequential stages:

  Stage 1  matching         - LCR and RCR both present (edit distance <= k),
                              read or its reverse complement.
  Stage 2  flanking         - exact left/right markers, extract what lies between.
  Stage 3  random region    - exact left/right stems with exactly N bases between.

Every stage (plus the raw input as "stage 0") is aggregated into a ranked table
(sequence<TAB>count, most frequent first, ties kept in first-seen order).

Stage 1 approximate matching uses edlib (compiled C/C++, Myers' bit-vector
algorithm) in infix ("HW") mode: the whole pattern must align, the matched
stretch of the read may start and end anywhere, substitutions / insertions /
deletions each cost 1 - the same semi-global edit distance as the original
pure-Python search, so the same reads pass.

Stage 1 prefilter (pigeonhole principle): each constant region is split once,
at start-up, into k+1 non-overlapping pieces that together cover it. k edits
can damage at most k of those pieces, so any occurrence within k edits leaves
at least one piece intact, i.e. present exactly in the read. If no piece
occurs exactly (fast substring search), the approximate match cannot succeed
and edlib is skipped. This only removes calls that would have returned False,
so the set of kept reads is unchanged.

Parallelism: reads are grouped into batches of --batch-size reads and the
per-read stage 1-3 work is done by a pool of worker processes. Workers only
compute; all counting happens in the main process, which consumes the batch
results strictly in submission order (= original read order). Every counter
therefore sees its sequences in exactly the same order as the single-process
version, so first-seen order - and with it the stable-sort tie order - is
unchanged.

Progress logging goes to stderr only (the Snakemake rule sends it to
logs/aptaselect.log); the result files are not affected by it.
"""

import argparse
import collections
import gzip
import io
import math
import multiprocessing as mp
import os
import sys
import time

import edlib

ORDINARY = frozenset("ACGT")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------------------------------------------------------------- logging
def log(msg):
    print(f"[aptaselect {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def fmt_int(n):
    return f"{n:,}"


def fmt_bytes(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n:.1f} TB"


def fmt_dur(sec):
    sec = int(round(max(0.0, sec)))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m{s:02d}s"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def pct(part, whole):
    return f"{100.0 * part / whole:.1f}%" if whole else "-"


def rss_bytes():
    """Resident memory of this (main) process, from /proc; None if unavailable."""
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        pass
    return None


def cgroup_memory_limit():
    """Container memory limit in bytes (cgroup v2 or v1), or None if unlimited."""
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(path) as fh:
                value = fh.read().strip()
            if value == "max":
                return None
            value = int(value)
            if value >= 1 << 60:  # v1 "unlimited" sentinel
                return None
            return value
        except (OSError, ValueError):
            continue
    return None


# --------------------------------------------------------------------------- IO
def open_fastq(path):
    """Open a FASTQ file for text reading, transparently handling gzip (detected
    by magic bytes). Returns (raw_binary_handle, text_handle); raw.tell() gives
    how far into the file on disk (compressed bytes for .gz) reading has got."""
    raw = open(path, "rb")
    magic = raw.read(2)
    raw.seek(0)
    if magic == b"\x1f\x8b":
        text = io.TextIOWrapper(gzip.GzipFile(fileobj=raw, mode="rb"))
    else:
        text = io.TextIOWrapper(raw)
    return raw, text


def iter_fastq_seqs(fh, path):
    """Yield the sequence line of every FASTQ record (upper-cased)."""
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
def split_pieces(pattern, max_err):
    """Split `pattern` into max_err+1 contiguous, non-overlapping pieces that
    cover it completely (sizes differ by at most 1). Computed once per primer."""
    n = max_err + 1
    base, extra = divmod(len(pattern), n)
    pieces, pos = [], 0
    for i in range(n):
        size = base + (1 if i < extra else 0)
        pieces.append(pattern[pos:pos + size])
        pos += size
    return tuple(pieces)


# per-process prefilter statistics (edlib calls made / skipped by the prefilter)
_STATS = [0, 0]


def approx_contains(pattern, pieces, text, max_err):
    """True if the whole `pattern` occurs somewhere in `text` with edit distance
    (substitution, insertion or deletion, each costing 1) <= max_err.

    `pieces` is split_pieces(pattern, max_err). If none of them occurs exactly
    in `text`, no occurrence within max_err edits can exist (pigeonhole), so
    return False without running the aligner. Otherwise edlib infix (HW)
    alignment decides: gaps before/after the matched stretch of text are free,
    the pattern must align end to end; with k=max_err edlib returns
    editDistance -1 when the best distance exceeds k.
    """
    m = len(pattern)
    if m == 0:
        return True
    if max_err >= m:
        return True  # the pattern can always be deleted entirely
    if max_err == 0:
        return pattern in text
    for piece in pieces:
        if piece in text:
            break
    else:
        _STATS[1] += 1
        return False  # no piece intact -> distance > max_err
    if not text:
        return False
    _STATS[0] += 1
    return edlib.align(pattern, text, mode="HW", task="distance", k=max_err)["editDistance"] != -1


def stage1_match(read, lcr, lcr_pieces, rcr, rcr_pieces, max_err):
    """Return the read (or its reverse complement) oriented LCR...RCR, else None."""
    if approx_contains(lcr, lcr_pieces, read, max_err) and approx_contains(rcr, rcr_pieces, read, max_err):
        return read
    rc = revcomp(read)
    if approx_contains(lcr, lcr_pieces, rc, max_err) and approx_contains(rcr, rcr_pieces, rc, max_err):
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


# ------------------------------------------------------------ per-batch worker
_P = None  # stage parameters, set once per worker process by _init_worker


def _init_worker(params):
    global _P
    _P = params


def process_batch(reads):
    """Run stages 1-3 on a list of reads; return (results, edlib_calls, skipped).

    results has one entry per read, same order: None (failed stage 1) or
    (oriented, flank, n_region) where flank / n_region are None if that stage
    failed. No counting of sequences here.
    """
    lcr, lcr_pieces, rcr, rcr_pieces, max_err, lm, rm, ls, rs, n_len = _P
    _STATS[0] = _STATS[1] = 0
    out = []
    append = out.append
    for read in reads:
        oriented = stage1_match(read, lcr, lcr_pieces, rcr, rcr_pieces, max_err)
        if oriented is None:
            append(None)
            continue
        flank = stage2_flank(oriented, lm, rm)
        if flank is None:
            append((oriented, None, None))
            continue
        append((oriented, flank, stage3_random_region(flank, ls, rs, n_len)))
    return out, _STATS[0], _STATS[1]


# ------------------------------------------------------------- CPU detection
def _cgroup_cpu_limit():
    """CPU limit from the cgroup CFS quota (containers), or None if unlimited."""
    # cgroup v2: cpu.max = "<quota|max> <period>"
    candidates = []
    try:
        with open("/proc/self/cgroup") as fh:
            for line in fh:
                parts = line.strip().split(":", 2)
                if len(parts) == 3 and parts[0] == "0" and parts[1] == "":
                    candidates.append(os.path.join("/sys/fs/cgroup", parts[2].lstrip("/"), "cpu.max"))
    except OSError:
        pass
    candidates.append("/sys/fs/cgroup/cpu.max")
    for path in candidates:
        try:
            with open(path) as fh:
                quota, period = fh.read().split()[:2]
            if quota == "max":
                return None
            return max(1, math.ceil(int(quota) / int(period)))
        except (OSError, ValueError):
            continue
    # cgroup v1: cpu.cfs_quota_us / cpu.cfs_period_us (quota -1 = unlimited)
    for base in ("/sys/fs/cgroup/cpu,cpuacct", "/sys/fs/cgroup/cpu"):
        try:
            with open(os.path.join(base, "cpu.cfs_quota_us")) as fh:
                quota = int(fh.read().strip())
            with open(os.path.join(base, "cpu.cfs_period_us")) as fh:
                period = int(fh.read().strip())
            if quota <= 0 or period <= 0:
                return None
            return max(1, math.ceil(quota / period))
        except (OSError, ValueError):
            continue
    return None


def usable_cpus():
    """CPUs actually available to this process: the CPU affinity mask (cpuset),
    further capped by any cgroup CPU quota. Not the machine total."""
    try:
        n = len(os.sched_getaffinity(0))
        src = "affinity"
    except (AttributeError, OSError):
        n = os.cpu_count() or 1
        src = "cpu_count"
    limit = _cgroup_cpu_limit()
    if limit is not None and limit < n:
        return limit, f"cgroup quota (affinity {n})"
    return max(1, n), src


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


def iter_batches(paths, batch_size):
    """Yield (file_index, reads, bytes_read_so_far_in_file) in file order then
    read order; each batch holds up to batch_size reads and never spans two files."""
    for idx, path in enumerate(paths):
        raw, fh = open_fastq(path)
        try:
            batch = []
            for read in iter_fastq_seqs(fh, path):
                batch.append(read)
                if len(batch) >= batch_size:
                    yield idx, batch, raw.tell()
                    batch = []
            if batch:
                yield idx, batch, raw.tell()
        finally:
            fh.close()
            raw.close()


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
    p.add_argument("--progress-every", type=int, default=100000,
                   help="log a progress line every this many reads (0 = off)")
    p.add_argument("--progress-seconds", type=float, default=60,
                   help="also log a progress line if this many seconds passed since the last one (0 = off)")
    p.add_argument("--workers", type=int, default=0,
                   help="worker processes; 0 = auto (CPUs allocated to this process, capped by --max-workers)")
    p.add_argument("--max-workers", type=int, default=0,
                   help="upper bound for auto-detected workers (e.g. Snakemake's thread allocation); 0 = no cap")
    p.add_argument("--batch-size", type=int, default=2000,
                   help="reads per work unit sent to a worker")
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
    if a.workers < 0:
        sys.exit("ERROR: workers must be >= 0 (0 = auto)")
    if a.batch_size < 1:
        sys.exit("ERROR: batch_size must be >= 1")
    for path in a.reads:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: input file not found: {path}")

    t_start = time.time()

    if a.workers > 0:
        workers, why = a.workers, "set explicitly"
    else:
        workers, why = usable_cpus()
        why = f"auto: {why}"
        if a.max_workers > 0 and a.max_workers < workers:
            workers, why = a.max_workers, f"{why}, capped by max_workers"

    # Primer piece split: computed once, shared by every read in every worker.
    lcr_pieces = split_pieces(lcr, a.max_errors)
    rcr_pieces = split_pieces(rcr, a.max_errors)
    insert_len = len(lcr) + len(left_stem) + a.random_length + len(right_stem) + len(rcr)

    sizes = [os.path.getsize(path) for path in a.reads]
    total_bytes = sum(sizes)
    mem_limit = cgroup_memory_limit()

    log("AptaSelect starting")
    log(f"settings: LCR={lcr} RCR={rcr} max_errors={a.max_errors} "
        f"markers={left_marker}/{right_marker} stems={left_stem}/{right_stem} "
        f"N={a.random_length} (insert length {insert_len})")
    log(f"prefilter pieces: LCR={lcr_pieces} RCR={rcr_pieces}")
    log(f"workers={workers} ({why}), batch_size={a.batch_size}, "
        f"memory limit={'none' if mem_limit is None else fmt_bytes(mem_limit)}")
    for i, (path, size) in enumerate(zip(a.reads, sizes), 1):
        log(f"input {i}/{len(a.reads)}: {path} ({fmt_bytes(size)})")
    log(f"output folder: {a.outdir}")

    params = (lcr, lcr_pieces, rcr, rcr_pieces, a.max_errors,
              left_marker, right_marker, left_stem, right_stem, a.random_length)

    os.makedirs(a.outdir, exist_ok=True)
    raw, s1, s2, s3 = OrderedCounter(), OrderedCounter(), OrderedCounter(), OrderedCounter()

    t0 = time.time()
    processed = 0
    stats = [0, 0]  # edlib calls, calls skipped by prefilter
    # progress state (touched only in consume, i.e. in read order)
    st = {
        "file": -1, "file_start_t": t0, "file_start_reads": 0,
        "done_bytes": 0,             # sizes of files already finished
        "last_t": t0, "last_reads": 0,
    }

    def progress_line(file_pos):
        now = time.time()
        elapsed = now - t0
        frac = (st["done_bytes"] + file_pos) / total_bytes if total_bytes else 0.0
        rate = processed / elapsed if elapsed > 0 else 0.0
        dt = now - st["last_t"]
        recent = (processed - st["last_reads"]) / dt if dt > 0 else 0.0
        eta = (elapsed / frac - elapsed) if frac > 0.001 else None
        rss = rss_bytes()
        log(f"{fmt_int(processed)} reads | {frac * 100:.1f}% of input | {fmt_dur(elapsed)} elapsed"
            f" | {fmt_int(int(rate))} reads/s (recent {fmt_int(int(recent))})"
            f" | ETA {fmt_dur(eta) if eta is not None else '-'}"
            f" | stage1 {fmt_int(s1.total)} ({pct(s1.total, processed)})"
            f" stage2 {fmt_int(s2.total)} ({pct(s2.total, processed)})"
            f" stage3 {fmt_int(s3.total)} ({pct(s3.total, processed)})"
            f" | unique raw {fmt_int(len(raw.counts))}"
            f" | main RSS {fmt_bytes(rss) if rss is not None else '-'}")
        st["last_t"], st["last_reads"] = now, processed

    def finish_file():
        idx = st["file"]
        if idx < 0:
            return
        n = processed - st["file_start_reads"]
        dt = time.time() - st["file_start_t"]
        log(f"finished file {idx + 1}/{len(a.reads)}: {fmt_int(n)} reads in {fmt_dur(dt)}")
        st["done_bytes"] += sizes[idx]

    def consume(file_idx, batch, file_pos, packed):
        """Count one batch in read order - the only place counters are touched."""
        nonlocal processed
        if file_idx != st["file"]:
            finish_file()
            st["file"] = file_idx
            st["file_start_t"] = time.time()
            st["file_start_reads"] = processed
            log(f"processing file {file_idx + 1}/{len(a.reads)}: {a.reads[file_idx]}")
        results, n_calls, n_skipped = packed
        stats[0] += n_calls
        stats[1] += n_skipped
        logged = False
        for read, r in zip(batch, results):
            processed += 1
            raw.add(read)
            if r is not None:
                oriented, flank, n_region = r
                s1.add(oriented)
                if flank is not None:
                    s2.add(flank)
                    if n_region is not None:
                        s3.add(n_region)
            if a.progress_every and processed % a.progress_every == 0:
                progress_line(file_pos)
                logged = True
        if not logged and a.progress_seconds and time.time() - st["last_t"] >= a.progress_seconds:
            progress_line(file_pos)

    if workers == 1:
        _init_worker(params)
        for file_idx, batch, pos in iter_batches(a.reads, a.batch_size):
            consume(file_idx, batch, pos, process_batch(batch))
    else:
        # Bounded, ordered dispatch: submit batches asynchronously, but always
        # collect the OLDEST pending batch first, so results are merged in the
        # original read order no matter which worker finishes first.
        max_inflight = workers * 4
        ctx = mp.get_context("fork")
        with ctx.Pool(workers, initializer=_init_worker, initargs=(params,)) as pool:
            pending = collections.deque()
            for file_idx, batch, pos in iter_batches(a.reads, a.batch_size):
                pending.append((file_idx, batch, pos, pool.apply_async(process_batch, (batch,))))
                while len(pending) >= max_inflight:
                    f, b, ps, res = pending.popleft()
                    consume(f, b, ps, res.get())
            while pending:
                f, b, ps, res = pending.popleft()
                consume(f, b, ps, res.get())
    finish_file()
    t_counted = time.time()
    log(f"all reads processed: {fmt_int(processed)} reads in {fmt_dur(t_counted - t0)}")
    checks = stats[0] + stats[1]
    if checks:
        log(f"stage1 approximate checks: {fmt_int(checks)}, edlib run: {fmt_int(stats[0])}, "
            f"skipped by prefilter: {fmt_int(stats[1])} ({pct(stats[1], checks)})")

    outputs = [
        (raw, "stage0_raw.ranked.tsv"),
        (s1, "stage1_matching.ranked.tsv"),
        (s2, "stage2_flanking.ranked.tsv"),
        (s3, "stage3_random_region.ranked.tsv"),
    ]
    for counter, name in outputs:
        tw = time.time()
        log(f"sorting and writing {name} ({fmt_int(len(counter.counts))} unique sequences)")
        write_ranked(counter, os.path.join(a.outdir, name))
        log(f"wrote {name} in {fmt_dur(time.time() - tw)}")

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
    log("wrote survival_summary.tsv")

    for name, value in summary:
        log(f"{name}\t{value}")
    t_end = time.time()
    rss = rss_bytes()
    log(f"timing: reading+processing {t_counted - t0:.1f}s, sorting+writing {t_end - t_counted:.1f}s, "
        f"total {t_end - t_start:.1f}s; main RSS at end {fmt_bytes(rss) if rss is not None else '-'}")
    log(f"done in {t_end - t0:.1f}s")


if __name__ == "__main__":
    main()