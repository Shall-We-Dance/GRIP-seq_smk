# Changelog

## 2.1.0

- Package the validated residual-adapter cleanup in a pinned Cutadapt environment; remove local executable and Python paths.
- Expose minimum retained read length and minimum aligned R2 length; optionally preserve raw concatenations without changing merge-before-dedup order.
- Add optional Input-enriched, DRACH-supported candidate m6A calling with explicit offsets, independent RT thresholds and an unselected positional diagnostic.
- Add generic library/group overlap tables, UpSet figures, directional containment, native-assembly Atlas BED6 validation and site-density metagenes.
- Report insufficient replication and zero-depth enrichment as not estimable; keep explicit QC exclusions in group definitions.
- Use equal-width UTR/CDS/UTR display in the example, and document signal CPM versus probability-density and candidate-site normalization.
- Retain vectorized deepTools 1-bp coverage output with numerical regression checks.
- Include synthetic algorithm/DAG tests, CI and a source-only privacy-checked release exporter.

Defaults in the new example are explicit; existing user configs are not overwritten.
The m6A module is optional and does not change motif-discovery thresholds.

## 2.0.1

- Correct the published GRIP library orientation: RNA is same as R2. Make this explicit and configurable; record it in signal statistics and provenance.
- Infer crosslinks one genomic base upstream of the R2 5-prime endpoint in RNA direction. Document unobserved upstream splice-junction limitations.
- Retain the crosslink shift margin in streaming counters to preserve coordinate ordering and count conservation.
- Add both-strand motif-coordinate, reference-boundary and adjacent-read ordering regression tests.
- Document the DRACH A +3/+4 convention and the distinct denominators in the published Fig. 6c comparison.

Critical correction: v2.0.0 strand-dependent tracks, crosslinks, motif conclusions and sense metaplots are invalid for the published GRIP protocol and require reanalysis. Whole-read BAM coverage is unaffected.

## 2.0.0

- Spliced, strand-aware 5UTR/CDS/3UTR metaplots with median-length display widths, probability density, per-gene signal and cohort QC.
- CIGAR-aware precise R2 endpoints and inferred crosslinks, RNA-strand count BEDs and tracks, clipping/blacklist audit.
- Unbiased RT site selection, reproducible bounded MEME pools and motif-offset/background statistics.
- Explicit BAM reanalysis input, library/reproducibility QC and private runtime provenance.
- Sparse interval BigWig arithmetic and synthetic algorithm regression tests.
- Source-only, history-free release export.

Breaking changes: site BED strands are RNA rather than R2; strict 5-prime endpoint filtering is the default; crosslinks are distinct from motif anchors; metaplot output names and profile schemas change.


## Earlier workflow updates

- Separate the public configuration template from the ignored local config.
- Add the Snakemake runner environment, synthetic DAG tests and GitHub Actions.
- Include STAR short-read settings and mapping-summary reporting.
- Include RT-termination motif discovery and TSS–TES metaplot modules.
- Explicitly select the complete workflow as the default target.
- Respect disabled IP/Input comparisons even when groups remain configured.
- Replace recursive temporary-directory cleanup with declared temporary outputs.
- Clarify sequence deduplication, R2 endpoints, strand handling and motif windows.

Existing installations retain their config.yaml and results. New clones must
copy config.example.yaml to config.yaml and supply reference/data paths before
running. Sequence deduplication is not UMI-aware. Analysis thresholds must be reported
from the effective configuration for each run.
