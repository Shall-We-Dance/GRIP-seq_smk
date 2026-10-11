# Analysis definitions (version 2.1.0)

## Relation to the published analysis

The GRIP-seq paper (Sun et al., Nature Chemistry 2023,
[doi:10.1038/s41557-022-01038-4](https://doi.org/10.1038/s41557-022-01038-4))
uses metaPlotR, MEME and R2-derived termination sites. The
[original metaplot script](https://github.com/Shall-We-Dance/GRIP-seq/blob/main/metaPlot/metaPlotR.r)
normalizes the three transcript segments and rescales UTR widths relative to
CDS. Version 2 follows this coordinate model and labels the differences below.

The publication mapped to hg19 and restricted analysis to CLIPper clusters.
This workflow accepts the supplied assembly and analyzes specified tracks or
count-defined RT local maxima; it does not currently execute CLIPper. The same
coordinate presentation does not imply an identical biological site set.
The original transcript-selection code overwrote its sorted table; this version
implements the documented longest-transcript intent deterministically.
Display medians use the common selected annotation cohort once, not a separate
sample-dependent, observation-weighted set. Profiles use exact binned signal
mass rather than an undocumented KDE bandwidth.

## Coordinates and 1-bp signal

All genomic positions are 0-based; BED intervals are half-open. In the published
GRIP library, R2 has the same orientation as RNA. The default is
`signal.rna_strand: same_as_read2`; every site BED denotes **RNA strand**.
Validate this relation using independent gene annotation for a new protocol.
Version 2.0.0 incorrectly assumed the opposite relation: its RNA-stranded tracks,
inferred crosslinks, motif results and sense metaplots must be recomputed.
Unstranded whole-read BAM coverage is unaffected by this correction.

The endpoint is the first query base aligned by CIGAR M/= /X from the sequencing
5′ end. Terminal D/N are not query bases. 5′ S/H/I means the physical RT boundary
is uncertain; the default rejects these reads. Terminal 3′ clipping alone is
allowed. Secondary, supplementary, unmapped, QC-failed and low-MAPQ alignments
are excluded. Proper-pair and duplicate flags have explicit optional policies.

The inferred crosslink is `signal.crosslink_shift` genomic bases upstream of
the R2 endpoint in RNA direction (default 1): forward R2 at p maps to p−1,
reverse R2 at p maps to p+1. This lies outside the observed R2 5′ boundary.
An upstream splice junction cannot be inferred from that R2 CIGAR; a genomic
adjacent coordinate at an exon boundary can therefore be intronic. This method
does not claim transcript-junction reconstruction. Out-of-reference coordinates
are rejected. Counts are buffered behind the alignment-start watermark by the
shift distance to preserve output order and aggregate nearby reads correctly.

An explicit `opposite_read2` mode supports other, independently validated
library orientations. Only in that mode does RNA-upstream advance into the
observed R2 CIGAR, skipping N gaps and counting reference D bases. The coordinate
mapping table and provenance always record the selected convention.

Each accepted R2 contributes one count. CPM = count / accepted R2 × 1,000,000;
both RNA strands share this denominator. The blacklist variant excludes an
observation if either its endpoint or inferred crosslink overlaps a merged
blacklist interval, and recomputes the denominator. Whole-read `bamCoverage`
uses its own aligned-read/blacklist normalization and is not numerically the
same quantity as R2 endpoint CPM. BigWigs with zero accepted R2 are valid empty
tracks accompanied by an explicit QC flag.

A crosslink is **not** a confirmed DRACH adenosine or m6A call. An optional
motif-conditioned diagnostic is disabled by default; it cannot feed discovery.
It does not shift already shifted coordinates a second time and reports
ambiguous equidistant motif matches.

## Spliced metaplot

Exons are merged and ordered in RNA 5′→3′ direction, removing introns. CDS plus
annotated stop-codon intervals define coding sequence. Models with missing
CDS, missing UTR, inconsistent coding boundaries or insufficient region length
are reported and excluded. Before choosing representatives, candidates must lie
on contigs shared by every configured metaplot sample/track and fit the shared
contig lengths; disagreeing lengths fail explicitly. Among these eligible
transcripts, choose longest spliced length, then longest CDS, then deterministic
transcript/locus order per gene. If a preferred patch-contig transcript is
unavailable, a compatible alternative is selected before excluding the gene.
The selection table reports these fallbacks and genes with no compatible
candidate. Only genes with all three regions enter the fixed shared cohort;
all samples use its common display widths.

Each region is divided into configurable relative bins (default 100). Signal is
integrated exactly over spliced genomic intervals and fractional bin boundaries.
Display width is median region length / median CDS length, or explicit fixed
widths. Probability density is total signal mass in each bin divided by total
mass and displayed bin width; its integral is 1 for nonzero profiles. This is
read-count weighted for endpoint/crosslink CPM. Mean CPM averages transcripts;
the additional equal-gene profile first divides each nonzero gene by its own
mean mature-transcript signal. Zero-signal genes have explicit counts.

RNA plus/minus tracks select sense signal. The BAM CPM panel is unstranded.
Overlapping genes can count the same genomic observation more than once;
per-gene tables and cohort selection make this visible. Replicate SD in grouped
plots is descriptive; one Input library supplies no biological-replicate SD.

## RT calling and motif positions

RT local maxima are called from unconditioned endpoint counts, requiring counts
strictly above `min_reads` (default 40) and at least `neighbor_fold` (default 1.5)
times both same-strand genomic neighbors. This is an explicit operational
interpretation of the published local-jump criterion, not a CLIPper replacement.
Splice-boundary neighborhood comparisons and genome-extracted sequence windows
are not transcriptome junction sequence reconstruction.

Discovery and validation are split by deterministic 10-kb genomic blocks shared
across strands, with 100-nt boundary guards; validation backgrounds stay within
the same block. Windows crossing a partition boundary are excluded.

MEME uses bounded, reproducibly selected loci, RNA-oriented windows and no
reverse-complement search. IP and Input pools are separate; repeated genomic
loci across libraries are deduplicated for pooled discovery. Configuration and
statistics record centering, discovery selection and holdout use.

For predefined motifs, an anchor is an explicitly specified motif position
(DRACH: middle A). Positive offset means downstream along RNA from the inferred
crosslink to the anchor. All matches are retained at each offset: one site can
match at several offsets, so sums across offsets need not be 1. Both independent
site frequencies and descriptive read-weighted frequencies are emitted.

Foreground windows are paired to valid nearby shifted genomic windows using a
fixed seed. These controls describe local genomic sequence background; they are
not expression-, GC-, or transcript-matched controls. Reported effective pair
counts exclude unusable windows. Paired discordant motif-presence counts are
compared with an exact one-sided binomial (McNemar) test, with BH correction
across tested motif/offset combinations and both crosslink/R2 populations per
sample. Site clustering and window
correlation limit inference; q values are exploratory positional enrichments,
not site-level modification probabilities. Input samples are analyzed separately
for comparison. MEME consensus scans are explicitly exploratory approximations
to motif PWMs, not FIMO PWM likelihood tests.

For comparison with published Fig. 6c, DRACH uses its central A as the motif
anchor. The paper places crosslinks at −3 (80.4%) and −4 (9.3%) relative to A;
our crosslink-centered axis therefore predicts A at +3 and +4. Those published
percentages are conditional on DRACH-containing reads, unlike our fraction of
all eligible RT sites at each offset. The denominators must not be conflated.

Peak offsets are estimated and reported, including ties and insufficient-data
cases, without changing crosslink coordinates or automatically assigning m6A.
Genome-derived sequence offsets describe contiguous genomic sequence in RNA
orientation, not a reconstructed transcript junction.

## Reanalysis and reproducibility

BAM input reruns filtering and all requested downstream modules. An optional
explicit imported BAM-coverage baseline may be reused only after record/filter/
normalization equivalence is checked; the input file identity is recorded. It does not
recreate reads discarded upstream and must not be used to claim validation of
changed preprocessing. FASTQ mode is required to re-evaluate trimming, sequence
deduplication, UMIs or alignment. The supplied sequence-dedup policy is not a
UMI-aware method. Runtime provenance includes effective config, input identities
and code hashes; real sample metadata never belongs in a public release.

## Optional candidate m6A inference and comparisons

`m6a` is separate from motif discovery and defaults to disabled. It calls RT
peaks independently with `count > min_reads` (default 10) and `count >= 1.5 *
max(left_neighbor, right_neighbor)`. The shared CIGAR-aware mapping gives each
peak's inferred crosslink. More than the configured minimum support must remain
after this mapping. No motif, database, local-background eligibility or MEME
partition is used to select RT peaks.

Input enrichment is `((C_IP+p)/N_IP) / ((C_Input+p)/N_Input)`, default `p=0.5` and
minimum fold 2. Counts come from the chosen blacklist variant; N is the total
accepted R2 **before blacklist filtering** in each library. This differs from
the denominator of a blacklist-specific CPM BigWig. Zero IP or Input depth is
reported as not estimable and yields no calls. This is an effect-size criterion,
not a modification probability, differential-binding test or FDR estimate.

The default adenosine is crosslink+3 in RNA-oriented genomic coordinates, and
must be the central A of DRACH (`[AGT][AG]AC[ACT]` in DNA notation). Disabling
`require_DRACH` retains only the A requirement and must be labeled accordingly.
Junction-spanning reference motifs are not reconstructed. Calls are
motif-supported inferred m6A candidates, not chemical confirmation. All RT
records and their exclusion flags remain available. The independent DRACH
position profile counts motifs within +/-25 nt of all sequence-eligible RT
sites **before** Input and motif selection, using crosslink-minus-A on the x axis.

Experimental sites are unique chromosome/adenosine-position/RNA-strand triples.
Group unions require >=1 included library; replicate sets require the configured
minimum support (default 2). Explicit exclusions are recorded. Insufficient
library counts are not estimable, not zero. UpSet columns are disjoint exact
membership patterns. Pairwise +/-1 and +/-5 analyses use maximum-cardinality
one-to-one matching. Directional containment is intersection/row set; Jaccard is
intersection/union. Empty-set fractions are NA.

Optional Atlas input must already be native-build BED6 single-adenosine sites.
All strand-oriented reference bases are validated and the file hash is recorded.
Exact matches and mutually exclusive distance categories (0, 1, 2-5, unmatched
within 5 nt) are reported. Database proximity is any-hit, unlike experimental
one-to-one matching. Main Atlas UpSets exclude Atlas-only background explicitly;
full-universe pattern tables are also exported. Database absence does not prove
novelty. Neither Atlas nor transcript assignment is required for a candidate call.

Candidate-site metagenes use representative complete spliced coding transcripts
selected from the configured GTF and available FASTA contigs. Each unique site
has total weight one, split over compatible genes. Each 5′ UTR/CDS/3′ UTR region
has equal display width with 30 bins by default. Density = weighted bin count /
(total assigned sites * scaled bin width), integrating to 1. Unassigned sites
remain in call catalogs. Signal metaplots remain a separate 100-bin calculation
from CPM tracks, with probability-density, mean-CPM and equal-gene panels.
