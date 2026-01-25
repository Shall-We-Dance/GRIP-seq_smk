# GRIP-seq Snakemake Analysis Pipeline

This repository provides a reproducible Snakemake workflow for **GRIP-seq** paired-end sequencing data. The pipeline performs per-lane QC with fastp, generates MultiQC summaries, aligns reads to a reference genome with STAR while retaining only uniquely mapped reads, and produces genome-browser-ready tracks including a **GRIP-seq 5′ end signal** derived from **Read 2 first base** (1 bp, strand-agnostic) with CPM normalization. BigWig tracks are generated both **with** and **without** blacklist filtering.

## Overview

### Inputs
- Paired-end raw FASTQ files (`*.fastq.gz`) for each sample.
- Multiple FASTQ pairs (lanes / technical replicates) may be assigned to the same sample.

### Outputs (per sample)
- **QC**
  - `results/qc/fastp/<sample>/unit<k>.html/json` (per-lane)
  - `results/qc/fastp/<sample>/fastp.html/json` (per-sample; computed after merge)
  - `results/qc/multiqc/multiqc_report.html`
- **Alignment**
  - `results/star/<sample>/<sample>.unique.mapq11.sorted.bam`
  - `results/star/<sample>/<sample>.unique.mapq11.sorted.bam.bai`
- **Signal tracks (BigWig, CPM normalized)**
  - BAM coverage:
    - `results/bigwig/<sample>/<sample>.bamCPM.noblacklist.bw`
    - `results/bigwig/<sample>/<sample>.bamCPM.blacklist.bw`
  - GRIP-seq 5′ end signal (Read2 first base, 1 bp):
    - `results/bigwig/<sample>/<sample>.R2firstbaseCPM.noblacklist.bw`
    - `results/bigwig/<sample>/<sample>.R2firstbaseCPM.blacklist.bw`
- **Genome browser track file**
  - `results/tracks/<sample>/<sample>.tracks.txt`

### Intermediate file handling
To minimize storage footprint, intermediate FASTQs produced by fastp and merged FASTQs are marked as temporary and are removed automatically by Snakemake. Final deliverables include:
- QC reports (fastp + multiqc)
- Final BAM + BAI
- BigWig tracks (coverage and 5′ end signal; with/without blacklist)
- Track text file

## Pipeline steps

1. **fastp QC (per-lane)**  
   Each input FASTQ pair is processed by fastp. Per-lane HTML/JSON reports are kept; cleaned FASTQs are temporary.

2. **Merge cleaned FASTQs**  
   Cleaned R1 files are concatenated to a single sample-level R1; cleaned R2 files are concatenated to a single sample-level R2 (gzip stream concatenation is valid).

3. **fastp (per-sample “report-only” pass)**  
   The merged FASTQs are passed through fastp with trimming/filtering disabled to produce a **sample-level report** consistent with the data used for alignment. Output FASTQs from this step are also temporary.

4. **STAR alignment (unique mapping)**  
   Reads are aligned with STAR using parameters that restrict to unique alignments (e.g., `--outFilterMultimapNmax 1`).

5. **Post-alignment MAPQ filtering**  
   The STAR BAM is further filtered with `MAPQ > 10` (implemented as `samtools view -q 11`). The resulting BAM is coordinate-sorted and indexed.

6. **BAM coverage BigWig (CPM)**  
   BigWigs are generated from the filtered BAM using deepTools `bamCoverage` with `--normalizeUsing CPM`, producing tracks with and without blacklist filtering.

7. **GRIP-seq 5′ end signal extraction (Read2 first base)**  
   The GRIP-seq 5′ signal is defined as the **first sequenced base of Read2**, projected onto the reference genome as a **1 bp** position per Read2 alignment:
   - If Read2 aligns to the **forward** strand: position = `reference_start`
   - If Read2 aligns to the **reverse** strand: position = `reference_end - 1`

   The signal is **strand-agnostic** (no +/- split) and is normalized to **CPM**, where the denominator is the number of usable mapped Read2 records after filtering (MAPQ and alignment flags). Tracks are generated with and without blacklist filtering.

8. **Blacklist filtering**  
   ENCODE blacklist BED files are downloaded from the Boyle-Lab blacklist repository and cached under `resources/blacklist/`. BigWigs are produced in paired versions (with/without blacklist) for both coverage and 5′ signal.

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

* `genome`: one of `mm10, mm39, hg19, hg38, ce10, ce11, dm3, dm6`
* `reference.star_index`: STAR genome index directory
* `reference.fasta`: reference FASTA (used to generate `.fai` / chrom sizes)
* `samples`: mapping of sample name to lists of FASTQs for R1 and R2

Example:

```yaml
genome: dm6

reference:
  star_index: "/path/to/STAR/index"
  fasta: "/path/to/genome.fa"
  gtf: "/path/to/genes.gtf"

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
* A Bowtie2 index path can also be specified in config for compatibility with broader projects, even though this workflow uses STAR.

## Running the workflow

From `workflow/`:

```bash
snakemake --use-conda --cores 16
```

To run a subset (example: one sample):

```bash
snakemake --use-conda --cores 16 results/bigwig/sampleA/sampleA.R2firstbaseCPM.blacklist.bw
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
* **Blacklist download issues**: Ensure network access to GitHub raw URLs; cached files are stored under `resources/blacklist/`.

## Citation / attribution

Blacklist BED files are retrieved from the Boyle-Lab blacklist repository (ENCODE blacklist lists). Please cite the appropriate sources in publications as required by the blacklist resource and your analysis conventions.


## Contact

For questions or extensions (e.g., track hubs, fragment-level scaling, additional QC metrics), please open an issue or contact the pipeline maintainer.
