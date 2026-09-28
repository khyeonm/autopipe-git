# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end SELEX FASTQ files. It runs three sequential filtering stages, then counts and ranks the sequences that survive each stage.

Insert layout: `LCR - left stem - N (random region) - right stem - RCR`

## Stages
Read 1 and Read 2 are each processed as their own sequence (all R1 reads first, then all R2 reads).

| Stage | What is kept |
|---|---|
| 0 raw | Every input read, as-is |
| 1 matching | Reads where the full LCR and the full RCR are each found anywhere in the read, within `lcr_max_errors` / `rcr_max_errors` edits (substitution, insertion or deletion). Their order is not checked. The read is checked as-is first, then as its reverse complement, and the orientation that carries both regions is kept. Otherwise the read is discarded. |
| 2 flanking | The sequence between the first `left_marker` and the farthest `right_marker` that can be reached through an uninterrupted run of A/C/G/T. Both markers must match exactly. There is no length limit and no stem check. If an N lies between the markers, or either marker is missing, the read yields nothing. |
| 3 random region | The region between exact `left_stem` and `right_stem` matches whose length is exactly `random_region_length`. Every overlapping left/right position pair is tried, and the first pair that satisfies the length wins. |

In each stage, identical sequences are counted and then sorted by count, highest first. The sort is stable, so sequences with equal counts stay in the order they were first seen.

## Inputs (config.yaml)
Every value is blank by default. You fill them in on the AutoPipe Input page.

| Key | Meaning |
|---|---|
| `r1`, `r2` | Paired-end FASTQ files (`.fastq` or `.fastq.gz`) |
| `lcr_seq`, `rcr_seq` | Left and right constant regions |
| `lcr_max_errors`, `rcr_max_errors` | Maximum edit distance allowed for each constant region |
| `left_marker`, `right_marker` | Stage 2 exact markers (at the inner end of the LCR and the inner start of the RCR) |
| `left_stem`, `right_stem` | Stage 3 stem sequences |
| `random_region_length` | Exact length of N |
| `threads` | Number of worker processes |

## Outputs (written directly to the output folder)
- `stage0_raw.ranked.tsv`, `stage1_matching.ranked.tsv`, `stage2_flanking.ranked.tsv`, `stage3_random_region.ranked.tsv`: each file has one header line (`sequence<TAB>count`), then one `sequence<TAB>count` row per unique sequence, sorted by count with the highest first.
- `survival_summary.tsv`: `name<TAB>value` rows with no header line. The rows are `total_reads_processed`, `raw_unique_sequences`, `stage1_matching_survivors`, `stage2_flanking_survivors`, `stage3_random_region_survivors` and `derived_insert_length`. The survivor values are total sequence counts, not unique counts. `derived_insert_length` = len(LCR) + len(left stem) + N + len(right stem) + len(RCR).
- `logs/aptaselect.log`

## Run
```bash
docker build -t aptaselect .
docker run --rm \
  -v /path/to/fastq:/input:ro \
  -v /path/to/results:/output \
  -v $(pwd)/config.yaml:/pipeline/config.yaml \
  aptaselect snakemake --cores 4
```