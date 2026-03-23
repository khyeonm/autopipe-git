#!/usr/bin/env python3
"""
AptaSelect: High-frequency aptamer candidate identification from paired-end FASTQ files.
Pipeline stages: Join -> Selection Filtering -> 1st Sort Filtering -> 2nd Sort Filtering -> Aggregation & Ranking
"""

import argparse
import gzip
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

# ─── Utilities ────────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def hamming(a, b):
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def phred_char_to_score(c):
    return ord(c) - 33

# ─── Stage 1: Join ────────────────────────────────────────────────────────────

def join_reads(r1_seq, r1_qual, r2_seq, r2_qual,
               min_overlap, pct_diff, is_short):
    rc_r2_seq  = reverse_complement(r2_seq)
    rc_r2_qual = r2_qual[::-1]

    if is_short:
        rd1_seq, rd1_qual = rc_r2_seq, rc_r2_qual
        rd2_seq, rd2_qual = r1_seq,    r1_qual
    else:
        rd1_seq, rd1_qual = r1_seq,    r1_qual
        rd2_seq, rd2_qual = rc_r2_seq, rc_r2_qual

    best_score = None
    best_ov    = None

    max_ov = min(len(rd1_seq), len(rd2_seq))
    for ov in range(min_overlap, max_ov + 1):
        tail = rd1_seq[-ov:]
        head = rd2_seq[:ov]
        d = hamming(tail, head)
        if d / ov <= pct_diff / 100.0:
            score = (1000 * (d * d + 1)) / ov
            if best_score is None or score < best_score:
                best_score = score
                best_ov    = ov

    if best_ov is None:
        return None, None

    ov = best_ov
    overlap_seq  = list(rd1_seq[-ov:])
    overlap_qual = list(rd1_qual[-ov:])
    for i in range(ov):
        q1 = phred_char_to_score(rd1_qual[-(ov - i)])
        q2 = phred_char_to_score(rd2_qual[i])
        if rd1_seq[-(ov - i)] != rd2_seq[i]:
            if q2 > q1:
                overlap_seq[i]  = rd2_seq[i]
                overlap_qual[i] = rd2_qual[i]

    joined_seq  = rd1_seq[:-ov] + "".join(overlap_seq) + rd2_seq[ov:]
    joined_qual = rd1_qual[:-ov] + "".join(overlap_qual) + rd2_qual[ov:]
    return joined_seq, joined_qual

# ─── Stages 2–4: Pattern Filtering ───────────────────────────────────────────

def match_pattern(seq, pattern, max_mismatches, start=0):
    plen = len(pattern)
    slen = len(seq)
    for pos in range(start, slen - plen + 1):
        if hamming(seq[pos:pos + plen], pattern) <= max_mismatches:
            yield pos

def filter_by_patterns(seq, pat1, pat2, max_mismatches, required_between=None):
    p1len = len(pat1)
    p2len = len(pat2)
    for pos1 in match_pattern(seq, pat1, max_mismatches):
        min_pos2 = pos1 + p1len
        for pos2 in match_pattern(seq, pat2, max_mismatches, start=min_pos2):
            between = pos2 - (pos1 + p1len)
            if required_between is None or between == required_between:
                return seq[pos1: pos2 + p2len]
    return None

# ─── FASTQ reading ────────────────────────────────────────────────────────────

def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def fastq_records(path):
    with open_fastq(path) as fh:
        while True:
            header = fh.readline().rstrip()
            if not header:
                break
            seq  = fh.readline().rstrip()
            _    = fh.readline()
            qual = fh.readline().rstrip()
            yield header, seq, qual

# ─── Chunk processing ─────────────────────────────────────────────────────────

def process_chunk(args):
    (chunk, min_overlap, pct_diff, is_short,
     sel_r1, sel_r2,
     s1_r1, s1_r2, s1_len,
     s2_r1, s2_r2, s2_len,
     max_mm) = args

    joined_seqs = []
    sel_seqs    = []
    sort1_seqs  = []
    sort2_seqs  = []

    for r1_seq, r1_qual, r2_seq, r2_qual in chunk:
        joined_seq, joined_qual = join_reads(
            r1_seq, r1_qual, r2_seq, r2_qual,
            min_overlap, pct_diff, is_short
        )
        if joined_seq is None:
            continue
        joined_seqs.append(joined_seq)

        sel = filter_by_patterns(joined_seq, sel_r1, sel_r2, max_mm, None)
        if sel is None:
            continue
        sel_seqs.append(sel)

        s1 = filter_by_patterns(sel, s1_r1, s1_r2, max_mm, s1_len)
        if s1 is None:
            continue
        sort1_seqs.append(s1)

        s2 = filter_by_patterns(s1, s2_r1, s2_r2, max_mm, s2_len)
        if s2 is None:
            continue
        sort2_seqs.append(s2)

    return (Counter(joined_seqs), Counter(sel_seqs),
            Counter(sort1_seqs),  Counter(sort2_seqs))

# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="AptaSelect pipeline")
    p.add_argument("--read1",  required=True)
    p.add_argument("--read2",  required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min_overlap",    type=int,   default=6)
    p.add_argument("--pct_diff",       type=float, default=8.0)
    p.add_argument("--is_short",       type=lambda x: x.lower() == "true", default=True)
    p.add_argument("--sel_read1",      default="CCACTTCTCCTTCCATCCTAAAC")
    p.add_argument("--sel_read2",      default="GAGTAGTTTGGAGGGTTGTCTG")
    p.add_argument("--s1_read1",       default="TCCTAAAC")
    p.add_argument("--s1_read2",       default="GAGTAGTT")
    p.add_argument("--s1_length",      type=int,   default=40)
    p.add_argument("--s2_read1",       default="TCTCTCTCTC")
    p.add_argument("--s2_read2",       default="GAGAGAGAGA")
    p.add_argument("--s2_length",      type=int,   default=20)
    p.add_argument("--max_mismatches", type=int,   default=1)
    p.add_argument("--top_n",          type=int,   default=10)  # kept for config compat, unused
    p.add_argument("--chunk_size",     type=int,   default=10000)
    p.add_argument("--threads",        type=int,   default=4)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pairs = list(zip(
        ((s, q) for _, s, q in fastq_records(args.read1)),
        ((s, q) for _, s, q in fastq_records(args.read2))
    ))
    all_records = [(r1s, r1q, r2s, r2q) for (r1s, r1q), (r2s, r2q) in pairs]

    chunk_size = args.chunk_size
    chunks = [all_records[i:i + chunk_size] for i in range(0, len(all_records), chunk_size)]

    task_args = [
        (chunk,
         args.min_overlap, args.pct_diff, args.is_short,
         args.sel_read1, args.sel_read2,
         args.s1_read1, args.s1_read2, args.s1_length,
         args.s2_read1, args.s2_read2, args.s2_length,
         args.max_mismatches)
        for chunk in chunks
    ]

    total_joined = Counter()
    total_sel    = Counter()
    total_sort1  = Counter()
    total_sort2  = Counter()

    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_chunk, a): i for i, a in enumerate(task_args)}
        for fut in as_completed(futures):
            jc, sc, s1c, s2c = fut.result()
            total_joined.update(jc)
            total_sel.update(sc)
            total_sort1.update(s1c)
            total_sort2.update(s2c)

    def write_tsv(counter, path, label):
        # All unique sequences, sorted by count descending, ranked from 1
        ranked = counter.most_common()
        with open(path, "w") as f:
            f.write("rank\tsequence\tcount\n")
            for rank, (seq, cnt) in enumerate(ranked, 1):
                f.write(f"{rank}\t{seq}\t{cnt}\n")
        print(f"[AptaSelect] {label}: {sum(counter.values())} total, "
              f"{len(counter)} unique -> all written to {path}", flush=True)

    write_tsv(total_joined, os.path.join(args.outdir, "joined_all.tsv"),    "Join")
    write_tsv(total_sel,    os.path.join(args.outdir, "selection_all.tsv"), "Selection")
    write_tsv(total_sort1,  os.path.join(args.outdir, "sort1_all.tsv"),     "1st Sort")
    write_tsv(total_sort2,  os.path.join(args.outdir, "sort2_all.tsv"),     "2nd Sort (final aptamers)")

    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write("AptaSelect Summary\n")
        f.write("==================\n")
        f.write(f"Total input pairs  : {len(all_records)}\n")
        f.write(f"Passed Join        : {sum(total_joined.values())}\n")
        f.write(f"Passed Selection   : {sum(total_sel.values())}\n")
        f.write(f"Passed 1st Sort    : {sum(total_sort1.values())}\n")
        f.write(f"Passed 2nd Sort    : {sum(total_sort2.values())}\n")

    print("[AptaSelect] Done.", flush=True)


if __name__ == "__main__":
    main()