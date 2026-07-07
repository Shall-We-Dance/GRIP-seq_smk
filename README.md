# GRIP-seq Snakemake Analysis Pipeline

This repository provides a reproducible Snakemake workflow for **GRIP-seq** paired-end sequencing data. The pipeline merges sample FASTQs, performs sample-level QC and trimming with fastp, generates MultiQC summaries, aligns reads to a reference genome with STAR while retaining only uniquely mapped reads, and produces genome-browser-ready tracks including a **GRIP-seq 5′ end signal** derived from **Read 2 first base** (1 bp, strand-agnostic) with CPM normalization. Optional blacklist-filtered BigWig tracks can be generated from a manually provided BED file or URL.

## Overview

### Inputs
- Paired-end raw FASTQ files (`*.fastq.gz`) for each sample.
- Multiple FASTQ pairs (lanes / technical replicates) may be assigned to the same sample.

### Outputs (per sample)
- **QC**
  - `results/qc/fastp/<sample>/merged_step1.html/json`
  - `results/qc/fastp/<sample>/merged_fastp_final.html/json`
  - `results/qc/multiqc/multiqc_report.html`
- **Alignment**
  - `results/star/<sample>/<sample>.unique.mapq11.sorted.bam`
  - `results/star/<sample>/<sample>.unique.mapq11.sorted.bam.bai`
- **Signal tracks (BigWig, CPM normalized)**
  - BAM coverage:
    - `results/bigwig/<sample>/<sample>.bamCPM.noblacklist.bw`
    - `results/bigwig/<sample>/<sample>.bamCPM.blacklist.bw` (if `filter_blacklist: true`)
  - GRIP-seq 5′ end signal (Read2 first base, 1 bp):
    - `results/bigwig/<sample>/<sample>.R2firstbaseCPM.noblacklist.bw`
    - `results/bigwig/<sample>/<sample>.R2firstbaseCPM.blacklist.bw` (if `filter_blacklist: true`)
- **Genome browser track file**
  - `results/tracks/<sample>/<sample>.tracks.txt`

### Intermediate file handling
To minimize storage footprint, intermediate FASTQs produced by fastp and merged FASTQs are marked as temporary and are removed automatically by Snakemake. Final deliverables include:
- QC reports (fastp + multiqc)
- Final BAM + BAI
- BigWig tracks (coverage and 5′ end signal; blacklist tracks are optional)
- Track text file

## Pipeline steps

1. **Merge raw FASTQs**
   R1 files are concatenated to a single sample-level R1; R2 files are concatenated to a single sample-level R2 (gzip stream concatenation is valid).

2. **fastp step 1: QC, adapter trimming, and optional dedup**
   The merged FASTQs are processed with fastp for QC, Illumina adapter trimming, and optional deduplication.

3. **fastp step 2: GRIP-seq custom trimming**
   The step 1 FASTQs are processed with fixed GRIP-seq trims. Output FASTQs are temporary and feed STAR alignment.

4. **STAR alignment (unique mapping)**  
   Reads are aligned with STAR using parameters that restrict to unique alignments (e.g., `--outFilterMultimapNmax 1`).

5. **Post-alignment MAPQ filtering**  
   The STAR BAM is further filtered with `MAPQ > 10` (implemented as `samtools view -q 11`). The resulting BAM is coordinate-sorted and indexed.

6. **BAM coverage BigWig (CPM)**
   BigWigs are generated from the filtered BAM using deepTools `bamCoverage` with `--normalizeUsing CPM`. If `filter_blacklist: true`, an additional blacklist-filtered track is generated.

7. **GRIP-seq 5′ end signal extraction (Read2 first base)**
   The GRIP-seq 5′ signal is defined as the **first sequenced base of Read2**, projected onto the reference genome as a **1 bp** position per Read2 alignment:
   - If Read2 aligns to the **forward** strand: position = `reference_start`
   - If Read2 aligns to the **reverse** strand: position = `reference_end - 1`

   The signal is **strand-agnostic** (no +/- split) and is normalized to **CPM**, where the denominator is the number of usable mapped Read2 records after filtering (MAPQ and alignment flags). If `filter_blacklist: true`, an additional blacklist-filtered track is generated.

8. **Blacklist filtering**
   Blacklist filtering is controlled by `filter_blacklist`. When enabled, provide exactly one of `blacklist.path` or `blacklist.url`. URL downloads are cached under `resources/blacklist/`.

## Requirements

- **Snakemake** (recommended: Snakemake >= 7)
- Conda / Mamba (recommended for environment management)
- STAR
- samtools
- fastp
- MultiQC
- deepTools
- Python packages: `pysam`, `pyBigWig`

All dependencies are provided via the conda environments under `workflow/envs/`.

## Installation

Recommended: create environments on-the-fly via Snakemake.

Example:
```bash
cd workflow
snakemake --use-conda --cores 16
```

To speed up conda solves, consider using mamba:

```bash
snakemake --use-conda --conda-frontend mamba --cores 16
```

## Configuration

Edit `workflow/config.yaml`.

Key fields:

* `genome`: genome label for reporting/logging
* `reference.star_index`: STAR genome index directory
* `reference.fasta`: reference FASTA (used to generate `.fai` / chrom sizes)
* `filter_blacklist`: enable/disable blacklist-filtered BigWig outputs
* `blacklist`: manual blacklist source; set exactly one of `path` or `url` when `filter_blacklist: true`
* `samples`: mapping of sample name to lists of FASTQs for R1 and R2

Example:

```yaml
genome: dm6

reference:
  star_index: "/path/to/STAR/index"
  fasta: "/path/to/genome.fa"
  gtf: "/path/to/genes.gtf"

filter_blacklist: true
blacklist:
  cache_dir: "resources/blacklist"
  # path: "/path/to/custom.blacklist.bed"
  # url: "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz"

samples:
  sampleA:
    R1:
      - "raw/sampleA_L001_R1.fastq.gz"
      - "raw/sampleA_L002_R1.fastq.gz"
    R2:
      - "raw/sampleA_L001_R2.fastq.gz"
      - "raw/sampleA_L002_R2.fastq.gz"
```

Notes:

* The workflow assumes paired-end reads and requires both R1 and R2 lists to be the same length per sample.

## Running the workflow

From `workflow/`:

```bash
snakemake --use-conda --cores 16
```

Dry run:

```bash
snakemake -n
```

## Outputs for genome browsers

BigWig files can be loaded directly into IGV or UCSC Genome Browser (if hosted).
A simple UCSC-style track file is generated at:

* `results/tracks/<sample>/<sample>.tracks.txt`

If using UCSC, update `bigDataUrl=` paths to valid HTTP(S) URLs pointing to your hosted BigWigs.

## Notes on normalization

* **BAM coverage tracks** use deepTools `bamCoverage --normalizeUsing CPM`.
* **GRIP-seq 5′ end tracks** are normalized by **CPM using Read2 counts** (one event per Read2), consistent with the Read2-first-base definition of the GRIP-seq 5′ signal.

## Troubleshooting

* **No reads after filtering**: If the 5′ extraction step fails with “No usable Read2 records”, confirm:

  * BAM contains properly paired alignments
  * Reads are mapped
  * MAPQ filtering is not overly stringent for your alignment settings
* **Blacklist configuration issues**: When `filter_blacklist: true`, set exactly one of `blacklist.path` or `blacklist.url`.
* **Blacklist download issues**: Ensure network access to the configured URL; cached files are stored under `resources/blacklist/`.

## Citation / attribution

If you use Boyle-Lab/ENCODE blacklist BED files, please cite the appropriate sources in publications as required by the blacklist resource and your analysis conventions.


## Contact

For questions or extensions (e.g., track hubs, fragment-level scaling, additional QC metrics), please open an issue or contact the pipeline maintainer.
