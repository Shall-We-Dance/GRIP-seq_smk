A paired-end GRIP-seq analysis workflow with fastp preprocessing, STAR alignment,
QC reports, Read 2 endpoint tracks, IP/Input comparisons, motif analyses, and
TSS–TES metaplots. The repository contains code and an example configuration;
sequencing data, experiment-specific settings, and analysis results stay local.

## Quick start

Run commands from the repository root on Linux, with Bash and Conda available.

```bash
git clone https://github.com/Shall-We-Dance/GRIP-seq_smk.git
cd GRIP-seq_smk
conda env create -f environment.yaml
conda activate gripseq-workflow

# New installations only: copy the template, then edit paths and sample groups.
cp -n config.example.yaml config.yaml

# Check the plan before running the analysis.
snakemake -s workflow/Snakefile --cores 16 --use-conda --dry-run
snakemake -s workflow/Snakefile --cores 16 --use-conda
```

The runner is pinned to Snakemake 9.3.3. Per-rule environments are specified in
`workflow/rules/envs/` and created by `--use-conda`. They pin the main tools,
but are not full platform-specific dependency lockfiles. Initial environment
creation and optional blacklist downloading require network access.

An existing local `config.yaml` remains usable after this repository update.
It is deliberately ignored by Git; publish `config.example.yaml` instead.
Relative data and reference paths are resolved from the working directory.

## Outputs

All paths below are relative to `output.dir`.

| Directory | Contents |
| --- | --- |
| `qc/fastp`, `qc/star`, `qc/multiqc` | Preprocessing and mapping QC in FASTQ mode |
| `star/<sample>` | Filtered coordinate-sorted BAM and BAI |
| `signal/<sample>` | RNA-strand endpoint/crosslink count BEDs, exact coordinate mapping and filtering statistics |
| `bigwig/<sample>` | Unstranded aligned-read CPM, plus combined and RNA plus/minus endpoint/crosslink CPM tracks |
| `bigwig/compare` | Group-average IP/Input comparisons |
| `metaplot` | Annotation selection, per-sample profiles, per-gene signal and grouped plots |
| `motif_de_novo/<sample>` | Sequence-unbiased RT peaks, sequence windows, MEME and selection statistics |
| `motif_de_novo/pooled/IP`, `pooled/Input` | Separate, locus-deduplicated discovery pools |
| `motif_distribution/<sample>` | Offset counts, local background comparisons, descriptive significance and figures |
| `motif_distribution/integrated.summary.tsv`, `grouped.predefined_motifs.pdf` | Per-sample motif offset summaries and IP/Input replicate comparisons |
| `qc/library_qc.*` | Library depth and precise-endpoint retention, plus gene-signal correlations |
| `provenance.json` | Private input/config/code manifest; BAM reanalysis inputs have SHA256 checksums |

`integrated.summary.tsv` includes every predefined and sample-specific de novo
motif, with its actual IUPAC consensus, anchor, effective denominator and q
values. Grouped curves compare predefined motifs only; same-named MEME models
with different consensuses are not pooled across samples.

## Key analysis choices

- **Precise endpoint filtering:** default `signal.five_prime_clip: exclude`
  rejects R2 5′ soft/hard clipping and insertion before the first aligned base.
  `aligned` retains the aligned endpoint for an explicitly labeled sensitivity
  analysis. Both modes report clipping; neither recovers a lost RT boundary.
- **Crosslinks:** one RNA nucleotide upstream of the R2 endpoint by default,
  following CIGAR splice gaps. Crosslinks, endpoints and putative m6A sites are
  distinct. Motif offsets are estimated without moving the crosslink coordinates.
- **Metaplots:** candidates are first restricted to compatible contigs shared
  across all configured sample/track BigWigs; then one longest eligible mature
  transcript is selected per gene, with explicit fallback/exclusion statistics.
  Exons are spliced, UTR/CDS/UTR regions are normalized separately, and display
  widths use the common cohort’s median segment lengths.
  Endpoint/crosslink tracks use the matching RNA strand. BAM coverage remains an
  unstranded comparison. Figures include signal probability density, mean CPM,
  and an additional equal-gene relative profile; endpoint/crosslink probability
  densities are read-count weighted. Group panels separate IP/Input.
- **Motifs:** RT peaks are selected by counts and neighbors before inspecting
  sequence. Predefined DRACH and held-out de novo consensus distributions are
  compared with nearby shifted genomic windows. The reports retain signed
  offsets, site and read-weighted frequencies, effective sample sizes, paired
  tests and BH-adjusted q values. These are exploratory positional statistics,
  not modification calls or biological-replicate differential-binding tests.
- **Sequence deduplication:** fastp `--dedup` is retained as an explicit option.
  Random library bases are trimmed, not extracted as UMIs; counts are not
  UMI-corrected molecule counts. Confirm library structure before changing trims.
- **QC:** low depth and excessive uncertain 5′ ends are flagged, without
  automatically removing samples. Correlations exclude shared-zero genes.

## References and license

GRIP/GECX-RNA method: Sun, W., Wang, N., Liu, H. et al.
[Genetically encoded chemical crosslinking of RNA in vivo](https://www.nature.com/articles/s41557-022-01038-4).
*Nature Chemistry* **15**, 21–32 (2023).

Related original analysis code: [Shall-We-Dance/GRIP-seq](https://github.com/Shall-We-Dance/GRIP-seq).
This workflow is an implementation with the choices described above, not a claim
of exact reproduction of every published analysis step.

MIT License — see [LICENSE](LICENSE).
