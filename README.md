# GRIP-seq Snakemake Analysis Pipeline

A reproducible Snakemake workflow for **GRIP-seq** (Genetically encoded RNA
Interactome Profiling / crosslinking + IP) paired-end sequencing data.

The pipeline merges per-sample FASTQs, performs two-step QC/trimming with
fastp, aligns with STAR (unique mappers only), filters alignments, and
produces genome-browser-ready tracks, including the GRIP-seq 5′ end signal
derived from the **Read 2 first base** (1 bp, CPM normalized). Optional
modules add blacklist-filtered tracks, deepTools `bigwigCompare`
log2(IP/Input) ratio tracks, and motif anchoring of the 1-nt crosslink sites
(e.g. m6A DRACH anchoring).

## Overview

### Inputs

- Paired-end raw FASTQ files (`*.fastq.gz`) per sample (from `config.yaml`).
- Multiple FASTQ pairs (lanes / technical replicates) may be assigned to the
  same sample.

### Outputs (per sample)

- **QC**
  - `results/qc/fastp/<sample>/merged_step1.html/json` (adapter trimming + dedup)
  - `results/qc/fastp/<sample>/merged_fastp_final.html/json` (GRIP-seq custom trimming)
  - `results/qc/multiqc/multiqc_report.html`
- **Alignment**
  - `results/star/<sample>/<sample>.unique.mapq11.sorted.bam` (+ `.bai`)
- **Signal tracks (BigWig, CPM normalized)**
  - BAM coverage: `results/bigwig/<sample>/<sample>.bamCPM[.blacklist].bw`
  - GRIP-seq 5′ end signal (Read2 first base, 1 bp): `results/bigwig/<sample>/<sample>.R2firstbaseCPM[.blacklist].bw`
- **Group comparisons (log2(IP/Input))**
  - `results/bigwig/compare/<group>.<track>.log2.bw` (per-sample BigWigs
    within each group side are averaged first, then compared)
- **Motif-anchored 1-nt sites**
  - `results/motif/<sample>/<sample>.anchored.bed` (BED6 + motif column)
- **Genome browser track file**
  - `results/tracks/<sample>/<sample>.tracks.txt`

### Intermediate file handling

Intermediate FASTQs and temporary files are marked as `temp()` and removed
automatically by Snakemake. Final deliverables are the QC reports, final
BAM + BAI, BigWig tracks, comparison tracks, anchored site BEDs and track
files.

## Pipeline steps

1. **Merge raw FASTQs** per sample (R1 and R2 concatenated separately).
2. **fastp step 1**: QC, Illumina adapter trimming, optional dedup
   (`fastp.dedup_adapter`).
3. **fastp step 2**: GRIP-seq custom fixed trimming (`fastp.grip_trim`).
   Disabling it makes step 2 a pass-through.
4. **STAR alignment**: unique alignments only
   (`--outFilterMultimapNmax 1`).
5. **Post-alignment filtering**: MAPQ > 10, explicit removal of secondary
   (0x100) and supplementary (0x800) alignments; coordinate-sorted + indexed.
6. **BAM coverage BigWig (CPM)** via deepTools `bamCoverage`; optional
   blacklist-filtered version.
7. **GRIP-seq 5′ end signal (Read2 first base)**: 1 bp position per Read2
   alignment (reference_start for forward, reference_end − 1 for reverse),
   strand-aware, CPM normalized by usable Read2 counts. An optional
   strand-aware 1-nt BED of these positions can be emitted.
8. **Blacklist filtering** (optional; URL-downloaded and cached under
   `resources/blacklist/`).
9. **bigwigCompare log2(IP/Input) tracks** (optional): per-group averaging of
   per-sample BigWigs, then `bigwigCompare --operation log2` with a
   configurable pseudocount and bin size.
10. **Motif anchoring** (optional): each R2-first-base 1-nt site is shifted to
    the cross-link site (Read2 is cDNA; `+site_shift` on '+' reads,
    `−site_shift` on '−' reads), optionally anchored to a motif on the RNA
    strand (default DRACH with the m6A at position 3), and a corrected 1-nt
    BED is written.

## Requirements

- Snakemake >= 7
- Conda / Mamba (recommended)
- Tools are installed on-the-fly from the pinned conda environments in
  `workflow/envs/` (`--use-conda`).

## Installation & usage

```bash
cd workflow
snakemake --use-conda --cores 16
```

Dry run:

```bash
snakemake -n --use-conda
```

To only run a specific module (e.g. motif anchoring):

```bash
snakemake --use-conda --cores 8 results/motif/sampleA_IP_1/sampleA_IP_1.anchored.bed
```

## Configuration

Edit `config.yaml`. Key fields:

- `reference.star_index` / `reference.fasta` / `reference.gtf`
- `filter_blacklist` and `blacklist` (exactly one of `path` / `url`)
- `samples`: sample name → lists of R1/R2 FASTQs
- `fastp.dedup_adapter` / `fastp.grip_trim` (enable flags respected)
- `star.multimap_nmax` (1 = unique alignments)
- `filtering.min_mapq` (default 11, i.e. MAPQ > 10)
- `bigwigCompare`: `enable`, `operation`, `pseudocount`, `binSize`,
  `tracks`, `groups` (IP/Input sample lists per comparison)
- `motif_anchoring`: `enable`, `motif`, `motif_pos`, `site_shift`,
  `search_window`, `keep_unmatched`

## Notes

- The workflow assumes paired-end reads; R1/R2 lists must have equal length
  per sample.
- The GRIP-seq 5′ end track is defined as the first sequenced base of Read2
  (1 bp per Read2), CPM normalized by Read2 counts.
- `bigwigCompare` only accepts two files, so replicates are averaged
  (per-base mean, missing = 0) before the log2 ratio is computed.
- At `binSize: 1` a whole-genome comparison takes ~1.5 h on hg38; use 10–25
  for fast browsing tracks.
- Motif anchoring is validated at workflow startup (motif length/position
  consistency, group sample names, etc.).

## License

MIT — see `LICENSE`.
