#!/usr/bin/env python3
"""
AptaSelect Stages 2-4: Sequential Pattern Validation.

Stage 2 (Selection): Find selection patterns, extract region (pattern-inclusive).
Stage 3 (1st Sort): Validate extracted sequence has sort patterns with correct between-length.
Stage 4 (2nd Sort): Validate extracted sequence has sort patterns with correct between-length.

Uses nested-loop matching: for each pattern1 match position, iterate all pattern2 
downstream positions. This avoids the greedy leftmost-only failure described in the spec.
"""

import sys
import json


def hamming_distance(s1, s2):
    """Compute Hamming distance between two equal-length strings."""
    return sum(1 for a, b in zip(s1, s2) if a != b)


def find_pattern_positions(seq, pattern, max_mismatches):
    """Find all positions where pattern matches seq with <= max_mismatches."""
    plen = len(pattern)
    positions = []
    for i in range(len(seq) - plen + 1):
        if hamming_distance(seq[i:i + plen], pattern) <= max_mismatches:
            positions.append(i)
    return positions


def find_pattern_pair(seq, pattern1, pattern2, max_mismatches, required_between=None):
    """
    Nested-loop search for a valid (pattern1, pattern2) pair.
    
    For each pattern1 match position, iterates all downstream pattern2 positions.
    If required_between is specified, checks that the region between patterns
    (excluding the patterns themselves) has exactly that length.
    
    Returns (p1_pos, p2_pos) or None if no valid pair found.
    """
    p1_len = len(pattern1)
    p2_len = len(pattern2)

    p1_positions = find_pattern_positions(seq, pattern1, max_mismatches)

    for p1_pos in p1_positions:
        # pattern2 must start after pattern1 ends
        search_start = p1_pos + p1_len
        for p2_pos in range(search_start, len(seq) - p2_len + 1):
            if hamming_distance(seq[p2_pos:p2_pos + p2_len], pattern2) <= max_mismatches:
                if required_between is not None:
                    between_length = p2_pos - (p1_pos + p1_len)
                    if between_length != required_between:
                        continue
                return (p1_pos, p2_pos)

    return None


def main():
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    config_path = sys.argv[3]

    with open(config_path) as f:
        config = json.load(f)

    sel_read1 = config['sel_read1']
    sel_read2 = config['sel_read2']
    s1_read1 = config['s1_read1']
    s1_read2 = config['s1_read2']
    s1_length = config['s1_length']
    s2_read1 = config['s2_read1']
    s2_read2 = config['s2_read2']
    s2_length = config['s2_length']
    max_mm = config['max_mismatches']

    total = 0
    passed_sel = 0
    passed_s1 = 0
    passed_s2 = 0

    with open(input_path) as inp, open(output_path, 'w') as out:
        for line in inp:
            seq = line.strip()
            if not seq:
                continue
            total += 1

            # Stage 2: Selection filtering
            result = find_pattern_pair(seq, sel_read1, sel_read2, max_mm, required_between=None)
            if result is None:
                continue
            p1_pos, p2_pos = result
            # Extract from start of pattern1 to end of pattern2 (inclusive)
            sel_seq = seq[p1_pos: p2_pos + len(sel_read2)]
            passed_sel += 1

            # Stage 3: 1st Sort filtering on sel_seq
            result = find_pattern_pair(sel_seq, s1_read1, s1_read2, max_mm, required_between=s1_length)
            if result is None:
                continue
            passed_s1 += 1

            # Stage 4: 2nd Sort filtering on sel_seq
            result = find_pattern_pair(sel_seq, s2_read1, s2_read2, max_mm, required_between=s2_length)
            if result is None:
                continue
            passed_s2 += 1

            out.write(sel_seq + '\n')

    # Write stats
    stats_path = output_path + '.stats'
    with open(stats_path, 'w') as sf:
        sf.write(f"input_sequences\t{total}\n")
        sf.write(f"passed_selection\t{passed_sel}\n")
        sf.write(f"passed_1st_sort\t{passed_s1}\n")
        sf.write(f"passed_2nd_sort\t{passed_s2}\n")

    print(f"Filtering: {total} input -> {passed_sel} selection -> {passed_s1} 1st sort -> {passed_s2} 2nd sort")


if __name__ == '__main__':
    main()