# GRIP-seq Snakemake workflow

**Correction:** version 2.0.0 assumed the wrong RNA/R2 strand relation. Recompute
its strand-dependent signals, crosslinks, motif analyses and metaplots with
version 2.1.0 (or the 2.0.1 correction); whole-read BAM coverage is unaffected.

Version 2.1.0 analyzes paired-end GRIP-seq libraries with explicit Read 2 endpoint
QC, RNA-strand crosslink tracks, segmented 5′UTR–CDS–3′UTR metaplots, and
motif-offset distributions with statistical reports. An optional module adds Input-enriched, motif-supported candidate m6A sites, replicate set comparisons, Atlas matching and site-density metagenes. Real sample configurations and results stay local.

## Run

```bash
git clone https://github.com/Shall-We-Dance/GRIP-seq_smk.git
cd GRIP-seq_smk
conda env create -f environment.yaml
conda activate gripseq-workflow
cp -n config.example.yaml config.yaml
# Edit references, samples and IP/Input groups first.
snakemake -s workflow/Snakefile --use-conda --cores 16 --resources mem_mb=64000 --dry-run
snakemake -s workflow/Snakefile --use-conda --cores 16 --resources mem_mb=64000 --rerun-incomplete
```

The runner is Snakemake 9.3.3. Rule environments pin their main tools, but are not
complete platform lockfiles. Conda environment creation and optional blacklist
downloads need network access. Use matching chromosome names and assembly for
FASTA, BAM, annotation and blacklist. Build the STAR index separately.

Each sample accepts matched nonempty `R1`/`R2` FASTQ lists (lanes from one
library), or a single `bam: "data/sample.bam"` for downstream reanalysis. Do not
supply both. BAM input is checked, filtered, sorted and indexed; original
FASTQ/STAR QC is not recreated. Previously removed reads cannot be recovered
from filtered BAMs. For changed trimming/alignment settings, use FASTQ input.
Biological replicates must remain separate samples.

In BAM mode only, `bam_coverage: {noblacklist: "tracks/sample.bw", blacklist:
"tracks/sample.filtered.bw"}` can import previously computed **aligned-read BAM
CPM coverage**. Use this only after verifying the same BAM records, 1-bp bins,
CPM normalization and blacklist; a matching chromosome header alone does not
establish that equivalence. Imported tracks are copied, and their hashes are
recorded in provenance. Endpoint/crosslink tracks are always computed anew.
When overriding an existing FASTQ config with `--configfile`, set `R1: []` and
`R2: []` in BAM samples to clear Snakemake's recursively merged defaults.

The GTF must contain transcript IDs, gene IDs, exons and CDS for the metaplot.
GTF `stop_codon` features are incorporated where present. Unsupported or
incomplete models are counted in the annotation report, not silently filled.
Reference FASTA needs an adjacent `.fai`, or a writable reference directory for
`samtools faidx`. Use a distinct output directory for a changed analysis and
preserve the prior configuration. Sample and group IDs use letters, digits,
underscores, hyphens and dots; avoid whitespace/shell metacharacters in paths.

## Outputs

All paths below are relative to `output.dir`.

| Directory | Contents |
| --- | --- |
| `qc/fastp`, `qc/cutadapt`, `qc/star`, `qc/multiqc` | Preprocessing and mapping QC in FASTQ mode |
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
| `m6a/<sample>` | All RT candidates, inferred m6A BED/TSV, Input enrichment, unselected DRACH offsets and calling parameters |
| `m6a/comparison` | Union/replicate sets, UpSet figures, containment, optional Atlas annotation and site metagenes |
| `provenance.json` | Private input/config/code manifest; BAM reanalysis inputs have SHA256 checksums |

See [METHODS.md](METHODS.md) for coordinates, statistical definitions and the
relationship to the published GRIP-seq analysis. All tabular data underlying
metaplots and motif distributions are retained.

`integrated.summary.tsv` includes every predefined and sample-specific de novo
motif, with its actual IUPAC consensus, anchor, effective denominator and q
values. Grouped curves compare predefined motifs only; same-named MEME models
with different consensuses are not pooled across samples.

## Key analysis choices

- **Precise endpoint filtering:** default `signal.five_prime_clip: exclude`
  rejects R2 5′ soft/hard clipping and insertion before the first aligned base.
  `aligned` retains the aligned endpoint for an explicitly labeled sensitivity
  analysis. Both modes report clipping; neither recovers a lost RT boundary.
- **Crosslinks:** R2 is RNA-sense by default (`signal.rna_strand: same_as_read2`).
  The crosslink is one genomic nucleotide upstream of its endpoint in RNA
  direction; an upstream splice junction outside the R2 alignment is unobserved. Crosslinks, endpoints and putative m6A sites are
  distinct. Motif offsets are estimated without moving the crosslink coordinates.
- **Metaplots:** candidates are first restricted to compatible contigs shared
  across all configured sample/track BigWigs; then one longest eligible mature
  transcript is selected per gene, with explicit fallback/exclusion statistics.
  Exons are spliced, UTR/CDS/UTR regions are normalized separately, and display
  widths are equal by default (`axis_scale: fixed`, 1:1:1); `median_length` is available for the metaPlotR-style display.
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

The historical BAM suffix `.unique.mapq11.sorted.bam` is retained for path
compatibility; the configured mapping/MAPQ settings are authoritative. Sparse
1-bp comparisons use the same per-base arithmetic with missing input treated as
zero, and explicitly retain zero output spans. Larger bins/other operations use
deepTools. Track declarations do not upload or host browser files.

## Candidate m6A analysis

Set `m6a.enable: true`, map each IP library to its Input in `m6a.controls`,
and list independent libraries in `m6a.groups`. This module has its own RT count
threshold (default strictly >10) and does not require MEME or metaplot modules.
Motif discovery retains its independently configured >40-read default. Compare
>20 and >40 in separate output directories when evaluating threshold sensitivity.

The default candidate requires a local RT peak, >=2-fold Input enrichment with
a 0.5-read pseudocount, and a DRACH central A at crosslink+3 along RNA. Each
library is called separately. These are inferred candidates; the Input fold is
an effect size, not a p value or FDR. DRACH is a selection condition here and must
not be presented as independent motif-enrichment evidence on the selected set.
Use the unselected-RT offset figure or the separate motif statistics for that.

The comparison exports group unions and sites supported by at least two
libraries. Groups with fewer eligible libraries are explicitly not estimable,
not empty biological sets. Review library QC first; optional `exclude_samples`
affects group comparisons only and is recorded in the output. There is no
implicit exclusion by sample name or automatic claim of biological replication.
The one-IP minimal template intentionally cannot estimate replicate support.

Optionally supply `m6a.atlas_bed` as BED6 with **0-based, half-open, 1-nt central
adenosines and RNA strand**, in the same genome build as the FASTA. Obtain and
convert the desired database release locally; its files are not bundled or
automatically downloaded. Every reference base is checked; no liftover is
performed. Set `atlas_label` to the release/build/snapshot used. Exact overlaps
are primary; +/-1 and +/-5 comparisons are also emitted. Experimental tolerance
matches are one-to-one; Atlas proximity is any-hit. Atlas-unmatched does not
establish a new modification. Main Atlas UpSets use the experimental universe;
full-universe pattern tables retain database-only sites explicitly.

Site metagenes give every unique candidate equal weight (fractional weights
for multiple compatible genes), use 30 bins per region, and integrate to one.
These differ from the 100-bin signal metaplots: only the Mean CPM signal panel
has CPM units; probability-density and equal-gene relative panels are unitless.
The CDS/3′ UTR boundary marks the end of the stop codon, not a fourth scaled region.

## Preprocessing portability

The example enables a pinned Cutadapt pass after fastp's GRIP-specific fixed
trims, with retained mates >=20 nt and R2 >=20 aligned M/=X bases. Adapter
sequences and thresholds are configurable. The residual-adapter audit samples
the first 20,000 retained pairs; it is not a random library-wide estimate.
`preprocessing.keep_merged_fastq: true` retains merged raw FASTQs for auditing;
the default removes these declared temporary copies after preprocessing to save
space. Original FASTQs are never changed. Merging always precedes deduplication.

BAM CPM coverage uses deepTools 3.5.6 with a vectorized 1-bp output writer; read
counting, filtering and normalization remain deepTools operations. Synthetic
regression tests compare its output with the unmodified bamCoverage executable.

## Validation and publication

```bash
# Add scientific test dependencies to the runner or use an analysis environment.
python -m pip install pysam pyBigWig numpy scipy matplotlib cutadapt==5.2 deepTools==3.5.6
python -m unittest discover -s tests -v
git diff --check
python tools/build_release.py --output dist/gripseq-v2.1.0
```

Tests exercise synthetic BAM CIGAR cases, RNA strand and coordinate conversion,
blacklist membership, exact BigWig math, spliced transcript models, motif-offset
statistics and Snakemake DAG variants. They complement real-data validation.

`tools/build_release.py` packages only explicitly selected code, documentation,
tests and the generic `config.example.yaml`; it omits local configs, `.git`
history, logs, data, references and results. The generated directory and archive
are suitable for review before copying into a clean Git checkout. Merely adding
`config.yaml` to `.gitignore` does not remove a previously tracked config or scrub
old Git history; inspect the destination repository before pushing.

## References

Sun W et al. [Genetically encoded chemical crosslinking of RNA in vivo](https://doi.org/10.1038/s41557-022-01038-4).
*Nature Chemistry* 15, 21–32 (2023).
[Original GRIP-seq analysis](https://github.com/Shall-We-Dance/GRIP-seq).
[metaPlotR](https://github.com/olarerin/metaPlotR).

MIT License; see [LICENSE](LICENSE).
