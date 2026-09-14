#!/usr/bin/env python3
"""
MEME motif analysis step for AptaSelect.

Takes the ranked random-region (N) table produced by AptaSelect
(stage3_random_region.ranked.tsv) — which is ALREADY trimmed to just the
variable random core, with primers, constant regions and stems removed — selects
the top-count cores, writes them to a FASTA file, and (optionally) runs MEME on
that FASTA to find the shared motif.

IMPORTANT: only the variable random core is fed to MEME. The constant/primer
regions are never included, because they would bias motif discovery (every
sequence shares them by construction).

Selection of "top" cores:
  --top-n N            : take the first N rows of the ranked table (default).
  --top-percent P      : if set (>0), take the top P% of the ranked list instead;
                         falls back to --top-n when unset / <= 0.
The ranked table is already sorted by count descending, so "top" = first rows.

This script always writes the selected cores to the FASTA output
(a real default output), regardless of whether MEME is actually run — so the
user can upload that FASTA to the MEME web server themselves.

MEME is run as the plain serial binary (no -p / MPI), called by full path,
treating sequences as DNA (-dna). Results go into <meme_out_dir>.

NOTE (not done here): a MEME motif is only a candidate. Afterwards the full
sequence should be checked separately for structure (G-quadruplex / stem-loop)
and folding energy to rule out motifs that are merely SELEX bias (e.g. cores
complementary to the constant regions). That structural check is a separate
step, not part of this pipeline.
"""

import argparse
import math
import os
import subprocess
import sys


def log(msg):
    print(f"[run_meme] {msg}", file=sys.stderr, flush=True)


def read_ranked_tsv(path):
    """Read AptaSelect ranked TSV (header: sequence<TAB>count), preserve order."""
    rows = []
    with open(path) as fh:
        header = fh.readline()  # skip 'sequence\tcount'
        if not header:
            return rows
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            seq = parts[0].strip()
            count = parts[1].strip() if len(parts) > 1 else ""
            if seq:
                rows.append((seq, count))
    return rows


def select_top(rows, top_n, top_percent):
    """Select the top rows. top_percent (>0) takes that share of the ranked
    list; otherwise fall back to top_n. The table is already count-sorted."""
    total = len(rows)
    if total == 0:
        return []
    if top_percent is not None and top_percent > 0:
        k = int(math.ceil(total * (top_percent / 100.0)))
        k = max(1, min(k, total))
        log(f"selecting top {top_percent}% -> {k} of {total} cores")
        return rows[:k]
    k = max(1, min(top_n, total))
    log(f"selecting top-n -> {k} of {total} cores")
    return rows[:k]


def write_fasta(rows, out_path):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as out:
        for i, (seq, count) in enumerate(rows, start=1):
            # header carries rank and count for traceability; MEME ignores it
            out.write(f">core{i}_count{count}\n{seq}\n")
    log(f"wrote {len(rows)} cores to {out_path}")


def run_meme(meme_bin, fasta_path, meme_out_dir, nmotifs, minw, maxw, mod,
             extra_args):
    """Run MEME as the plain serial binary. Never pass -p (MPI fails in the
    container)."""
    cmd = [
        meme_bin,
        fasta_path,
        "-dna",                 # nucleotide alphabet
        "-oc", meme_out_dir,    # overwrite/create output dir (keeps re-runs clean)
        "-nmotifs", str(nmotifs),
        "-mod", mod,
    ]
    if minw is not None and int(minw) > 0:
        cmd += ["-minw", str(minw)]
    if maxw is not None and int(maxw) > 0:
        cmd += ["-maxw", str(maxw)]
    if extra_args:
        cmd += extra_args
    log("running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    log(f"MEME finished; results in {meme_out_dir}")


def main():
    ap = argparse.ArgumentParser(description="MEME motif step for AptaSelect")
    ap.add_argument("--ranked-tsv", required=True,
                    help="AptaSelect stage3 random-region ranked TSV")
    ap.add_argument("--fasta-out", required=True,
                    help="Path to write selected top cores (always written)")
    ap.add_argument("--meme-out-dir", required=True,
                    help="Directory for MEME results (meme_out)")
    ap.add_argument("--top-n", type=int, default=100,
                    help="Number of top cores to select (fallback)")
    ap.add_argument("--top-percent", type=float, default=0.0,
                    help="If >0, take this %% of the ranked list instead of top-n")
    ap.add_argument("--run-meme", type=lambda s: str(s).lower() in
                    ("1", "true", "yes", "on"), default=False,
                    help="Actually run MEME (default off: only write FASTA)")
    ap.add_argument("--meme-bin", default="/opt/meme/bin/meme",
                    help="Full path to the serial MEME binary")
    ap.add_argument("--nmotifs", type=int, default=3)
    ap.add_argument("--minw", type=int, default=0)
    ap.add_argument("--maxw", type=int, default=0)
    ap.add_argument("--mod", default="zoops",
                    help="MEME site distribution: oops | zoops | anr")
    ap.add_argument("--meme-extra", default="",
                    help="Extra raw args passed through to MEME")
    args = ap.parse_args()

    rows = read_ranked_tsv(args.ranked_tsv)
    log(f"read {len(rows)} ranked cores from {args.ranked_tsv}")
    if not rows:
        log("WARNING: ranked table is empty; writing empty FASTA")

    selected = select_top(rows, args.top_n, args.top_percent)
    write_fasta(selected, args.fasta_out)

    if not args.run_meme:
        log("run_meme is OFF — FASTA written, MEME not run. "
            "Upload the FASTA to the MEME web server, or enable run_meme.")
        return

    if not selected:
        log("ERROR: no cores selected; cannot run MEME")
        sys.exit(1)

    if not os.path.isfile(args.meme_bin) or not os.access(args.meme_bin, os.X_OK):
        log(f"ERROR: MEME binary not found/executable at {args.meme_bin}")
        sys.exit(1)

    extra = args.meme_extra.split() if args.meme_extra.strip() else []
    run_meme(args.meme_bin, args.fasta_out, args.meme_out_dir,
             args.nmotifs, args.minw, args.maxw, args.mod, extra)


if __name__ == "__main__":
    main()