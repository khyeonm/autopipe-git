#!/usr/bin/env python3
"""
AptaSelect Stage 5: Aggregation & Ranking.
Counts identical sequences, sorts by frequency, outputs ALL sequences sorted descending.
Also outputs top N separately.
"""

import sys
import json
from collections import Counter


def main():
    input_path = sys.argv[1]
    output_all_path = sys.argv[2]
    output_top_path = sys.argv[3]
    config_path = sys.argv[4]

    with open(config_path) as f:
        config = json.load(f)

    top_n = config['top_n']

    # Count all sequences
    counts = Counter()
    with open(input_path) as inp:
        for line in inp:
            seq = line.strip()
            if seq:
                counts[seq] += 1

    # Sort by count descending (all sequences)
    ranked_all = counts.most_common()

    # Write ALL sequences sorted by count
    with open(output_all_path, 'w') as out:
        out.write("rank\tsequence\tcount\n")
        for rank, (seq, count) in enumerate(ranked_all, 1):
            out.write(f"{rank}\t{seq}\t{count}\n")

    # Write top N
    ranked_top = ranked_all[:top_n]
    with open(output_top_path, 'w') as out:
        out.write("rank\tsequence\tcount\n")
        for rank, (seq, count) in enumerate(ranked_top, 1):
            out.write(f"{rank}\t{seq}\t{count}\n")

    total_unique = len(counts)
    total_seqs = sum(counts.values())
    print(f"Aggregation: {total_seqs} sequences, {total_unique} unique, top {min(top_n, len(ranked_top))} reported")

    # Write stats
    stats_path = output_top_path + '.stats'
    with open(stats_path, 'w') as sf:
        sf.write(f"total_sequences\t{total_seqs}\n")
        sf.write(f"unique_sequences\t{total_unique}\n")
        sf.write(f"top_n_reported\t{min(top_n, len(ranked_top))}\n")


if __name__ == '__main__':
    main()