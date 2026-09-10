# GRIP-seq Snakemake workflow

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

## Inputs and configuration

Edit `config.yaml` before running:

- `reference.star_index`: a prebuilt STAR genome index compatible with the
  configured STAR version. Index building is not part of this workflow.
- `reference.fasta`: matching genome FASTA. Its adjacent `.fai` must exist,
  or the reference directory must be writable so `samtools faidx` can create it.
- `reference.gtf`: matching GTF/GFF annotation, required for metaplots.
- `samples`: non-empty R1/R2 lists of gzip-compressed FASTQs. Both lists must
  have the same length and corresponding file order. Multiple lanes from the
  same library may be combined; keep independent biological libraries separate.
- `bigwigCompare.groups`: explicit IP and Input sample lists for each comparison.
- `output.dir`: output directory (default `results`). Add custom output
  directories to `.gitignore` if they reside inside the repository.

Use consistent chromosome names across the index, FASTA, annotation and blacklist.
The example blacklist is for hg38; change it for other assemblies. Use paths
without whitespace or shell metacharacters with the current shell rules.

The example includes two placeholder samples, not downloadable demonstration
data. The small synthetic test inputs described below are only for DAG checks.

## Workflow and outputs

1. Merge FASTQ lanes per library, then run fastp adapter/QC processing and
   optional sequence deduplication.
2. Run a second fastp pass for library-specific fixed trimming.
3. Align with STAR, filter primary alignments by MAPQ, sort and index BAMs.
4. Generate QC summaries, alignment coverage and R2 endpoint CPM BigWigs.
5. Run enabled comparison, motif and metaplot modules.

| Output under `results/` | Contents |
| --- | --- |
| `qc/fastp/<sample>/` | Both fastp HTML/JSON reports |
| `qc/star/short_read_mapping_summary.tsv` | STAR mapping/retention metrics |
| `qc/multiqc/multiqc_report.html` | Aggregated QC |
| `star/<sample>/` | Filtered BAM/BAI, STAR log and splice-junction table |
| `bigwig/<sample>/` | `bamCPM.noblacklist.bw`, `R2firstbaseCPM.noblacklist.bw`, optional `.blacklist.bw` variants |
| `bigwig/compare/` | Comparisons after averaging each group side's CPM tracks |
| `motif/<sample>/` | Specified-motif anchored sites |
| `motif_de_novo/<sample>/` | RT-site BED, sequence FASTA, statistics and MEME outputs |
| `motif_de_novo/integrated/` | Pooled sequences and MEME outputs |
| `metaplot/` | Per-sample profiles and combined TSS–TES plots |
| `tracks/<sample>/` | Browser track declarations; hosting BigWigs is a separate step |

BAM filenames retain the historical `.unique.mapq11.sorted.bam` suffix even
if `filtering.min_mapq` is changed; the actual threshold is the configured value.
The same applies to the `unique` label if `star.multimap_nmax` is increased.
Snakemake manages declared `temp()` outputs. The workflow does not recursively
delete the output temporary directory; STAR may leave additional auxiliary files.

## Analysis definitions and limitations

### Trimming and duplicates

The example preserves the existing GRIP preprocessing settings: R1 front/tail
12/10 nt, R2 front/tail 10/12 nt, and fastp deduplication enabled. Review these
against the actual library structure, read length and adapter read-through.

`fastp.dedup_adapter.dedup: true` invokes fastp sequence deduplication.
This release does **not** extract UMIs or perform UMI-aware BAM deduplication.
Setting it to false disables that step without adding a replacement.
The fixed front trims discard potential barcode bases. For UMI analysis, preserve
those bases before trimming and add a separately validated UMI workflow.
See the [fastp duplicate definition](https://github.com/OpenGene/fastp#duplication-rate-and-deduplication).

`dedup_adapter.enable: false` disables deduplication, not the first QC pass.
`grip_trim.enable: false` disables fixed trims and adapter trimming in the
second pass, but fastp's normal quality/length filters still run.

### Alignment and R2 signal

The supplied STAR settings relax aligned-length/score fractions for short inserts
while retaining absolute match and mismatch thresholds. These are analysis
choices to review with QC, not a universal CLIP-seq preset.

The R2 endpoint is `reference_start` for forward alignments and
`reference_end - 1` for reverse alignments. This is the first **aligned** base
in sequencing orientation; 5′ soft-clipped bases are not recovered.
R2 CPM uses the count of usable R2 alignments after the selected filtering.
BigWigs sum both strands; the motif BED retains strand-specific counts.
The BAM CPM track is aligned-read coverage from `bamCoverage`, not a deduplicated
molecule count or explicitly reconstructed full-fragment coverage.

### Optional modules

- **Blacklist**: `filter_blacklist` controls the extra filtered tracks. Specify
  exactly one of `blacklist.path` or `blacklist.url` when enabled.
- **IP/Input**: set `bigwigCompare.enable`. Replicate CPM tracks are averaged
  with missing values treated as zero, then compared with a pseudocount.
  These are visualization tracks, not replicate-aware differential-binding tests.
- **Specified motifs**: `motif_anchoring.enable` shifts R2 endpoints and can
  select or reposition them relative to a supplied motif. This is separate from
  unbiased motif discovery; review `search_window` before interpretation.
- **De novo motifs**: `de_novo_motif.enable` calls strand-specific R2-start
  positions with counts strictly greater than `min_reads` and at least
  `neighbor_fold` times each adjacent position. It scans genome-wide and does
  not run CLIPper cluster calling. BED sites are shifted one nucleotide toward
  the inferred crosslink; FASTA windows remain centered on the unshifted R2 start
  and are oriented in the inferred RNA direction. MEME runs without reverse-
  complement searching. Low-sequence inputs generate an explicit “not run”
  report and placeholder logo. `samples: all` pools IP and Input alike;
  supply an IP-only list if that is the intended integrated analysis.
- **Metaplots**: `metaplot.enable` orients genes by annotation strand, scales
  gene bodies and plots both unstranded signal tracks. This does not separate
  sense/antisense RNA signal. `samples: all` includes every configured sample.

Module parameters are documented inline in `config.example.yaml`.
To request one output explicitly:

```bash
snakemake -s workflow/Snakefile --use-conda --cores 8 \
  results/motif/example_IP_1/example_IP_1.anchored.bed
```

## Development and validation

```bash
conda activate gripseq-workflow
python -m unittest discover -s tests -v
git diff --check
```

The tests parse Python/YAML and build synthetic Snakemake DAGs for all modules,
core-only mode, disabled comparisons, and disabled blacklist filtering.
They check the default target and shell-command expansion without executing
bioinformatics tools. GitHub Actions runs the same checks. They do not validate
biological accuracy, install every rule environment, or replace a real-data run.

Repository layout:

```text
config.example.yaml         Public configuration template
environment.yaml            Snakemake runner environment
workflow/Snakefile          Default all target and module includes
workflow/rules/             Processing rules
workflow/rules/scripts/     Python analysis helpers
workflow/rules/envs/        Per-rule Conda environments
tests/                      Offline synthetic DAG checks
.github/workflows/ci.yaml   Automated checks
```

## References and license

GRIP/GECX-RNA method: Sun, W., Wang, N., Liu, H. et al.
[Genetically encoded chemical crosslinking of RNA in vivo](https://www.nature.com/articles/s41557-022-01038-4).
*Nature Chemistry* **15**, 21–32 (2023).

Related original analysis code: [Shall-We-Dance/GRIPseq](https://github.com/Shall-We-Dance/GRIP-seq).
This workflow is an implementation with the choices described above, not a claim
of exact reproduction of every published analysis step.

MIT License — see [LICENSE](LICENSE).
