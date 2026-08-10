#!/usr/bin/env python3
"""
extract_cores.py -- Motif-analysis pre-step for AptaSelect.

Takes the top-ranked sequences from an AptaSelect ranked TSV, locates the two
experiment-specific flanking motifs, and writes ONLY the variable core (the
region strictly between the two flanks) to a FASTA file. Primers and other
fixed regions are excluded -- this FASTA is what MEME (or the MEME website)
consumes.

How many top-ranked sequences are used:
  * If top_percent is set (non-blank), use that top PERCENTAGE of the ranked
    sequences (highest-ranked first). This overrides top_n.
  * If top_percent is blank, fall back to top_n (a fixed count).

Hard requirement: the two flank sequences and the core length (gap) have NO
defaults. If any is missing/blank, the script exits non-zero and writes nothing.
"""

import argparse
import sys
import os
import math


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def find_flank(seq, flank, max_mm):
    """Return (start, end_index_just_after) of the first match of `flank` in
    `seq`, or (-1, -1). With max_mm=0 this is an exact search."""
    n, m = len(seq), len(flank)
    for i in range(n - m + 1):
        window = seq[i:i + m]
        if max_mm == 0:
            if window == flank:
                return i, i + m
        else:
            if hamming(window, flank) <= max_mm:
                return i, i + m
    return -1, -1


def read_all_sequences(tsv_path):
    """Read an AptaSelect ranked TSV (rank<TAB>count<TAB>sequence) and return
    ALL (rank, count, sequence) tuples in rank order (highest-ranked first)."""
    rows = []
    with open(tsv_path) as fh:
        fh.readline()  # header: rank  count  sequence
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rows.append((parts[0], parts[1], parts[2]))
    return rows


def main():
    p = argparse.ArgumentParser(description="Extract variable cores between two flanking motifs.")
    p.add_argument("--input", required=True, help="AptaSelect ranked TSV")
    p.add_argument("--output", required=True, help="Output FASTA of variable cores")
    p.add_argument("--flank5", default="", help="5' flanking motif (immediately before core)")
    p.add_argument("--flank3", default="", help="3' flanking motif (immediately after core)")
    p.add_argument("--gap", default="", help="Exact core length in bp between the two flanks")
    p.add_argument("--top-n", type=int, default=100,
                   help="Fixed number of top-ranked sequences to use (fallback when --top-percent is blank)")
    p.add_argument("--top-percent", default="",
                   help="If set, use this top percentage of the ranked sequences (highest-ranked first); overrides --top-n")
    p.add_argument("--flank-max-mismatches", type=int, default=0)
    args = p.parse_args()

    # -- Hard guard: no defaults. Refuse to run on blank motif/gap values. --
    flank5 = (args.flank5 or "").strip().upper()
    flank3 = (args.flank3 or "").strip().upper()
    gap_raw = (args.gap or "").strip()

    missing = []
    if not flank5:
        missing.append("core_flank5 (5' flanking motif)")
    if not flank3:
        missing.append("core_flank3 (3' flanking motif)")
    if gap_raw == "":
        missing.append("core_gap (variable core length)")
    if missing:
        sys.stderr.write(
            "ERROR: motif analysis requires values that have no defaults.\n"
            "Missing/blank: " + ", ".join(missing) + "\n"
            "Fill these in config.yaml (or pass them in) and re-run. Nothing was written.\n"
        )
        sys.exit(2)

    try:
        gap = int(gap_raw)
    except ValueError:
        sys.stderr.write(f"ERROR: core_gap must be an integer bp length, got '{gap_raw}'. Nothing was written.\n")
        sys.exit(2)
    if gap < 0:
        sys.stderr.write("ERROR: core_gap must be >= 0. Nothing was written.\n")
        sys.exit(2)

    if not os.path.exists(args.input):
        sys.stderr.write(f"ERROR: input table not found: {args.input}\n")
        sys.exit(1)

    all_rows = read_all_sequences(args.input)
    total = len(all_rows)

    # -- Decide how many top-ranked sequences to use --
    # top_percent (when set) overrides top_n; otherwise fall back to top_n.
    top_percent_raw = (args.top_percent or "").strip()
    if top_percent_raw != "":
        try:
            pct = float(top_percent_raw)
        except ValueError:
            sys.stderr.write(
                f"ERROR: top_percent must be a number (a percentage), got '{top_percent_raw}'. Nothing was written.\n"
            )
            sys.exit(2)
        if pct <= 0:
            sys.stderr.write(f"ERROR: top_percent must be greater than 0, got {pct:g}. Nothing was written.\n")
            sys.exit(2)
        if pct > 100:
            sys.stderr.write(f"NOTE: top_percent {pct:g} > 100; clamping to 100 (using all ranked sequences).\n")
            pct = 100.0
        k = math.ceil(total * pct / 100.0)
        if k > total:
            k = total
        selection_desc = f"top {pct:g}% of {total} ranked sequences"
    else:
        k = args.top_n
        selection_desc = f"top_n = {args.top_n} (top_percent blank)"

    selected = all_rows[:k]

    n_written = n_no5 = n_no3 = n_badlen = 0

    outdir = os.path.dirname(os.path.abspath(args.output))
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    with open(args.output, "w") as out:
        for rank, count, seq in selected:
            s = seq.strip().upper()
            _, l5_end = find_flank(s, flank5, args.flank_max_mismatches)
            if l5_end < 0:
                n_no5 += 1
                continue
            rel_start, _ = find_flank(s[l5_end:], flank3, args.flank_max_mismatches)
            if rel_start < 0:
                n_no3 += 1
                continue
            core = s[l5_end:l5_end + rel_start]
            if len(core) != gap:
                n_badlen += 1
                continue
            out.write(f">rank{rank}_count{count}_len{len(core)}\n{core}\n")
            n_written += 1

    sys.stderr.write(
        "=== extract_cores summary ===\n"
        f"Input table:              {args.input}\n"
        f"Ranked sequences (total): {total}\n"
        f"Selection:                {selection_desc}\n"
        f"Sequences selected:       {len(selected)}\n"
        f"5' flank:                 {flank5}\n"
        f"3' flank:                 {flank3}\n"
        f"Core length (gap):        {gap}\n"
        f"Cores written:            {n_written}\n"
        f"Dropped (no 5' flank):    {n_no5}\n"
        f"Dropped (no 3' flank):    {n_no3}\n"
        f"Dropped (len != gap):     {n_badlen}\n"
        f"Output FASTA:             {args.output}\n"
    )

    if n_written == 0:
        sys.stderr.write(
            "WARNING: no cores matched the given flanks/length. The FASTA is empty. "
            "Check core_flank5, core_flank3, and core_gap against your library design.\n"
        )


if __name__ == "__main__":
    main()