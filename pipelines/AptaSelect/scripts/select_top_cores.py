#!/usr/bin/env python3
"""Write the top-count random-region (N) cores as FASTA for MEME.

Input is AptaSelect's final ranked table (stage3_random_region.ranked.tsv):
a header line "sequence<TAB>count", then one row per unique N sequence, most
frequent first. Those sequences are already just the variable core - primers,
constant regions and stems were removed by AptaSelect - so they are written
as-is; nothing is added or trimmed here.

How many are taken:
  * --top-percent P (0 < P <= 100) given: ceil(P% of the unique sequences in
    the ranked table), at least 1;
  * otherwise: the first --top-n rows (or all rows if there are fewer).
The cut is by rank in the table's own order; if the last taken row has the
same count as the next row, that tie is split by the table's order and a
warning is logged.

FASTA header: >rank<R>_count<C> (rank 1 = most frequent).
If the table has no sequences, an empty FASTA is written and a warning logged,
so the count-sorting run still finishes; the MEME step refuses empty input.
"""

import argparse
import math
import sys

ORDINARY = frozenset("ACGT")


def log(msg):
    print(f"[select_top_cores] {msg}", file=sys.stderr, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True, help="stage3_random_region.ranked.tsv")
    p.add_argument("--out", required=True, help="output FASTA")
    p.add_argument("--top-n", type=int, required=True)
    p.add_argument("--top-percent", type=float, default=None)
    a = p.parse_args()

    if a.top_n < 1:
        sys.exit("ERROR: meme_top_n must be >= 1")
    if a.top_percent is not None and not (0 < a.top_percent <= 100):
        sys.exit("ERROR: meme_top_percent must be > 0 and <= 100 (leave it blank to use meme_top_n)")

    rows = []
    with open(a.table) as fh:
        header = fh.readline().rstrip("\n")
        if header != "sequence\tcount":
            sys.exit(f"ERROR: unexpected header in {a.table}: {header!r} (expected 'sequence<TAB>count')")
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                sys.exit(f"ERROR: {a.table} line {lineno}: expected 2 tab-separated columns")
            seq, count = parts[0].strip().upper(), int(parts[1])
            rows.append((seq, count))

    n = len(rows)
    if a.top_percent is not None:
        k = max(1, math.ceil(n * a.top_percent / 100.0)) if n else 0
        how = f"top {a.top_percent:g}% of {n:,} unique sequences"
    else:
        k = min(a.top_n, n)
        how = f"top {a.top_n:,} (of {n:,} unique sequences)"

    selected = rows[:k]
    written = 0
    skipped = 0
    with open(a.out, "w") as out:
        for rank, (seq, count) in enumerate(selected, start=1):
            if not seq or set(seq) - ORDINARY:
                skipped += 1
                log(f"WARNING: rank {rank} skipped (empty or non-ACGT): {seq!r}")
                continue
            out.write(f">rank{rank}_count{count}\n{seq}\n")
            written += 1

    log(f"selection: {how} -> {k:,} rows, {written:,} written to {a.out}"
        + (f", {skipped} skipped" if skipped else ""))
    if selected:
        lengths = sorted({len(s) for s, _ in selected})
        log(f"counts {selected[0][1]:,} .. {selected[-1][1]:,}; core length(s): {lengths}; "
            f"reads covered: {sum(c for _, c in selected):,} of {sum(c for _, c in rows):,}")
    if 0 < k < n and rows[k - 1][1] == rows[k][1]:
        tied = sum(1 for _, c in rows if c == rows[k - 1][1])
        log(f"WARNING: cut falls inside a tie - {tied} sequences share count {rows[k - 1][1]}; "
            f"kept by table order")
    if written == 0:
        log("WARNING: no sequences written - the random-region table is empty; MEME cannot run on this")


if __name__ == "__main__":
    main()