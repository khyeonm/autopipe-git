#!/usr/bin/env python3
"""AptaSelect: rank high-frequency aptamer candidates from paired-end SELEX FASTQ.

Stages (each read from R1 and R2 is handled as its own sequence):
  0  raw          every input read, as-is
  1  matching     LCR and RCR both present (edit distance <= tolerance, anywhere,
                  no order requirement); read kept as-is, else reverse complement
                  kept if it carries both; else discarded
  2  flanking     sequence between first left marker and the farthest right marker
                  reachable through an uninterrupted A/C/G/T run (exact markers)
  3  random N     region between left stem and right stem of exactly the target
                  length; all overlapping stem position pairs are tried, first valid
                  pair wins

Every stage's surviving sequences are counted and stably sorted by count only
(descending), so ties keep first-seen input order.
"""
import argparse
import gzip
import os
import re
import sys
from collections import Counter
from multiprocessing import Pool

try:
    import edlib  # fast Myers bit-vector edit distance
    HAVE_EDLIB = True
except ImportError:  # pragma: no cover - pure-Python fallback
    HAVE_EDLIB = False

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def _contains_py(pattern, text, k):
    """Sellers semi-global DP: does pattern occur anywhere in text with <= k edits
    (substitution, insertion or deletion)?"""
    m = len(pattern)
    if m == 0:
        return True
    prev = list(range(m + 1))  # column for empty text prefix
    if prev[m] <= k:
        return True
    for c in text:
        cur = [0] * (m + 1)  # free start anywhere in text
        for i in range(1, m + 1):
            cost = 0 if pattern[i - 1] == c else 1
            cur[i] = min(prev[i - 1] + cost, prev[i] + 1, cur[i - 1] + 1)
        if cur[m] <= k:
            return True
        prev = cur
    return False


def contains(pattern, text, k):
    if HAVE_EDLIB:
        # mode HW = infix: pattern must align fully, gaps at text ends are free
        d = edlib.align(pattern, text, mode="HW", task="distance", k=k)["editDistance"]
        return 0 <= d <= k  # edlib ignores k for empty text, so re-check the bound
    return _contains_py(pattern, text, k)


# ---------------------------------------------------------------- worker state
P = {}


def init_worker(params):
    global P
    P = dict(params)
    P["flank_re"] = re.compile(
        re.escape(P["left_marker"]) + "([ACGT]*)" + re.escape(P["right_marker"])
    )


def all_positions(s, sub):
    """All (overlapping) start positions of sub in s."""
    out, i = [], s.find(sub)
    while i != -1:
        out.append(i)
        i = s.find(sub, i + 1)
    return out


def stage1(read):
    lcr, rcr = P["lcr"], P["rcr"]
    if contains(lcr, read, P["lcr_err"]) and contains(rcr, read, P["rcr_err"]):
        return read
    rc = revcomp(read)
    if contains(lcr, rc, P["lcr_err"]) and contains(rcr, rc, P["rcr_err"]):
        return rc
    return None


def stage2(seq):
    start = seq.find(P["left_marker"])  # first (leftmost) left marker
    if start == -1:
        return None
    m = P["flank_re"].match(seq, start)  # greedy run -> farthest reachable right marker
    if m is None:
        return None
    return m.group(1)  # no length limit, no stem check


def stage3(seq):
    ls, rs, n = P["left_stem"], P["right_stem"], P["n_len"]
    lefts = all_positions(seq, ls)
    if not lefts:
        return None
    rights = set(all_positions(seq, rs))
    for i in lefts:  # all left/right pairs; the length constraint fixes the partner
        j = i + len(ls) + n
        if j in rights:
            return seq[i + len(ls): j]
    return None


def process_chunk(reads):
    out = []
    for r in reads:
        s1 = stage1(r)
        s2 = stage2(s1) if s1 is not None else None
        s3 = stage3(s2) if s2 is not None else None
        out.append((s1, s2, s3))
    return out


def process_chunk_pair(reads):
    return reads, process_chunk(reads)


# ---------------------------------------------------------------------- I/O
def read_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            yield seq.upper()


def chunked(paths, size):
    buf = []
    for p in paths:
        for seq in read_fastq(p):
            buf.append(seq)
            if len(buf) >= size:
                yield buf
                buf = []
    if buf:
        yield buf


def write_ranked(counter, path, header):
    # dict/Counter keeps first-seen order; sorted() is stable -> ties stay in that order
    ranked = sorted(counter.items(), key=lambda kv: kv[1], reverse=True)
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for seq, c in ranked:
            fh.write(f"{seq}\t{c}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", default="")
    ap.add_argument("--r2", default="")
    ap.add_argument("--lcr", default="")
    ap.add_argument("--rcr", default="")
    ap.add_argument("--lcr-max-errors", default="")
    ap.add_argument("--rcr-max-errors", default="")
    ap.add_argument("--left-marker", default="")
    ap.add_argument("--right-marker", default="")
    ap.add_argument("--left-stem", default="")
    ap.add_argument("--right-stem", default="")
    ap.add_argument("--random-region-length", default="")
    ap.add_argument("--header", default="sequence\tcount")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--chunk-size", type=int, default=20000)
    a = ap.parse_args()

    def clean(x):
        x = "" if x is None else str(x).strip()
        return "" if x.lower() in ("none", "null") else x

    required = {
        "r1": a.r1, "r2": a.r2, "lcr_seq": a.lcr, "rcr_seq": a.rcr,
        "lcr_max_errors": a.lcr_max_errors, "rcr_max_errors": a.rcr_max_errors,
        "left_marker": a.left_marker, "right_marker": a.right_marker,
        "left_stem": a.left_stem, "right_stem": a.right_stem,
        "random_region_length": a.random_region_length,
    }
    missing = [k for k, v in required.items() if clean(v) == ""]
    if missing:
        sys.exit("ERROR: these config.yaml values are blank and must be filled in: "
                 + ", ".join(missing))

    def as_int(name, v):
        try:
            iv = int(clean(v))
        except ValueError:
            sys.exit(f"ERROR: {name} must be a non-negative integer, got {v!r}")
        if iv < 0:
            sys.exit(f"ERROR: {name} must be a non-negative integer, got {v!r}")
        return iv

    params = {
        "lcr": clean(a.lcr).upper(), "rcr": clean(a.rcr).upper(),
        "lcr_err": as_int("lcr_max_errors", a.lcr_max_errors),
        "rcr_err": as_int("rcr_max_errors", a.rcr_max_errors),
        "left_marker": clean(a.left_marker).upper(),
        "right_marker": clean(a.right_marker).upper(),
        "left_stem": clean(a.left_stem).upper(),
        "right_stem": clean(a.right_stem).upper(),
        "n_len": as_int("random_region_length", a.random_region_length),
    }
    for f in (a.r1, a.r2):
        if not os.path.isfile(clean(f)):
            sys.exit(f"ERROR: input FASTQ not found: {f}")

    derived_insert_length = (len(params["lcr"]) + len(params["left_stem"]) + params["n_len"]
                             + len(params["right_stem"]) + len(params["rcr"]))

    raw, s1c, s2c, s3c = Counter(), Counter(), Counter(), Counter()
    n_raw = n1 = n2 = n3 = 0
    paths = [clean(a.r1), clean(a.r2)]  # all R1 reads, then all R2 reads

    def consume(chunk, results):
        nonlocal n_raw, n1, n2, n3
        for read, (s1, s2, s3) in zip(chunk, results):
            raw[read] += 1
            n_raw += 1
            if s1 is not None:
                s1c[s1] += 1
                n1 += 1
            if s2 is not None:
                s2c[s2] += 1
                n2 += 1
            if s3 is not None:
                s3c[s3] += 1
                n3 += 1

    print(f"edlib available: {HAVE_EDLIB}; threads: {a.threads}", file=sys.stderr)
    if a.threads > 1:
        with Pool(a.threads, initializer=init_worker, initargs=(params,)) as pool:
            # imap returns chunks in input order, so first-seen order is preserved
            for ch, res in pool.imap(process_chunk_pair, chunked(paths, a.chunk_size)):
                consume(ch, res)
    else:
        init_worker(params)
        for ch in chunked(paths, a.chunk_size):
            consume(ch, process_chunk(ch))

    os.makedirs(a.outdir, exist_ok=True)
    o = a.outdir
    write_ranked(raw, os.path.join(o, "stage0_raw.ranked.tsv"), a.header)
    write_ranked(s1c, os.path.join(o, "stage1_matching.ranked.tsv"), a.header)
    write_ranked(s2c, os.path.join(o, "stage2_flanking.ranked.tsv"), a.header)
    write_ranked(s3c, os.path.join(o, "stage3_random_region.ranked.tsv"), a.header)

    summary = [
        ("total_reads_processed", n_raw),
        ("raw_unique_sequences", len(raw)),
        ("stage1_matching_survivors", n1),
        ("stage2_flanking_survivors", n2),
        ("stage3_random_region_survivors", n3),
        ("derived_insert_length", derived_insert_length),
    ]
    with open(os.path.join(o, "survival_summary.tsv"), "w") as fh:
        for k, v in summary:
            fh.write(f"{k}\t{v}\n")
    for k, v in summary:
        print(f"{k}\t{v}", file=sys.stderr)


if __name__ == "__main__":
    main()