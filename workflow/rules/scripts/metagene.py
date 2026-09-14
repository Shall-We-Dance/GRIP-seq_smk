"""Spliced 5'UTR/CDS/3'UTR annotation and exact sparse BigWig integration.

Coordinates are zero-based, half-open. Segment blocks are stored in RNA 5'->3'
order; no intronic bases enter a region. The coding segment includes explicitly
annotated stop codons (GTF producers disagree whether CDS includes these bases).
"""

from collections import Counter, defaultdict
import gzip
import json
import re
import statistics

REGIONS = ("utr5", "cds", "utr3")
REGION_LABELS = ("5′UTR", "CDS", "3′UTR")
SELECTION_POLICY = (
    "Among transcripts with valid contiguous spliced CDS and all three regions "
    "at least min_region_length nt, select one per gene by longest spliced transcript, "
    "longest CDS including stop codon, then lexical transcript ID, chromosome "
    "and exon coordinates. Explicit cds_start_NF/cds_end_NF annotations are excluded."
)


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if start >= end:
            raise ValueError("Empty/reversed interval")
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def parse_attributes(text):
    """GTF attributes, including repeated tag fields (GFF3 is not silently guessed)."""
    attrs = defaultdict(list)
    for match in re.finditer(r'(?:^|;)\s*([^\s;]+)\s+"([^"]*)"', text):
        attrs[match.group(1)].append(match.group(2))
    return attrs


def _spliced_intervals(blocks, exons, strand):
    """Map genomic blocks into transcript coordinates; reject nonexonic bases."""
    mapped = []
    offset = 0
    covered = 0
    for start, end in (exons if strand == "+" else exons[::-1]):
        for a, b in blocks:
            left, right = max(start, a), min(end, b)
            if left >= right:
                continue
            covered += right - left
            mapped.append([offset + (left - start if strand == "+" else end - right),
                           offset + (right - start if strand == "+" else end - left)])
        offset += end - start
    if covered != sum(b - a for a, b in blocks):
        raise ValueError("coding_feature_outside_exons")
    return merge_intervals(mapped)


def _genomic_blocks(left, right, exons, strand):
    blocks = []
    offset = 0
    for start, end in (exons if strand == "+" else exons[::-1]):
        a, b = max(left, offset), min(right, offset + end - start)
        if a < b:
            if strand == "+":
                blocks.append([start + a - offset, start + b - offset])
            else:
                blocks.append([end - b + offset, end - a + offset])
        offset += end - start
    return blocks


def segment_transcript(record, min_region_length=1):
    """Validate a coding transcript and return its mature RNA region blocks."""
    if min_region_length < 1:
        raise ValueError("min_region_length must be >= 1")
    if record.get("partial_cds", False):
        raise ValueError("partial_cds_annotation")
    if not record.get("exons"):
        raise ValueError("no_exons")
    if not record.get("cds"):
        raise ValueError("noncoding")
    exons = merge_intervals(record["exons"])
    cds = merge_intervals(record["cds"])
    spliced_cds = _spliced_intervals(cds, exons, record["strand"])
    if len(spliced_cds) != 1:
        raise ValueError("discontinuous_spliced_cds")
    coding_start, coding_end = spliced_cds[0]
    stop = record.get("stop_codons", [])
    if stop:
        mapped_stop = _spliced_intervals(merge_intervals(stop), exons, record["strand"])
        # A split codon is contiguous after exon concatenation. A CDS-inclusive
        # stop overlaps the last three coding bases; a CDS-exclusive stop starts
        # exactly at coding_end. Reject unrelated internal/remote stop features.
        if (len(mapped_stop) != 1 or mapped_stop[0][1] - mapped_stop[0][0] != 3
                or mapped_stop[0][0] not in (coding_end - 3, coding_end)):
            raise ValueError("inconsistent_stop_codon")
        coding_end = max(coding_end, mapped_stop[0][1])
    total = sum(b - a for a, b in exons)
    lengths = [coding_start, coding_end - coding_start, total - coding_end]
    if any(length < min_region_length for length in lengths):
        raise ValueError("missing_or_short_utr_or_cds")
    edges = [0, coding_start, coding_end, total]
    result = {key: record[key] for key in ("gene_id", "transcript_id", "chrom", "strand")}
    result["gene_name"] = record.get("gene_name", record["gene_id"])
    result["lengths"] = dict(zip(REGIONS, lengths))
    result["regions"] = {
        region: _genomic_blocks(edges[i], edges[i + 1], exons, record["strand"])
        for i, region in enumerate(REGIONS)
    }
    result["exons"] = exons
    result["stop_codon_annotated"] = bool(stop)
    return result


def prepare_annotation(gtf_path, min_region_length=1, available_chrom_sizes=None):
    """Build one deterministic representative per gene and an auditable cohort."""
    records = {}
    parse_counts = Counter()
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF row {number}: expected 9 columns")
            chrom, _, feature, start, end, _, strand, _, text = fields
            if feature.lower() not in ("transcript", "exon", "cds", "stop_codon"):
                continue
            attrs = parse_attributes(text)
            gene_id = attrs.get("gene_id", [None])[0]
            tx_id = attrs.get("transcript_id", [None])[0]
            if not gene_id or not tx_id:
                parse_counts["rows_missing_gtf_gene_or_transcript_id"] += 1
                continue
            if strand not in ("+", "-"):
                parse_counts["rows_invalid_strand"] += 1
                continue
            a, b = int(start) - 1, int(end)
            if a < 0 or b <= a:
                raise ValueError(f"Invalid GTF coordinates at row {number}")
            # RefSeq can reuse transcript IDs on alternate contigs. They must
            # remain distinct until the gene representative is selected.
            key = (gene_id, tx_id, chrom, strand)
            record = records.setdefault(key, {
                "gene_id": gene_id, "transcript_id": tx_id, "chrom": chrom,
                "strand": strand, "gene_name": attrs.get("gene_name", [gene_id])[0],
                "exons": [], "cds": [], "stop_codons": [], "partial_cds": False,
            })
            record["partial_cds"] |= bool(
                {"cds_start_NF", "cds_end_NF"} & set(attrs.get("tag", []))
            )
            destination = {"exon": "exons", "cds": "cds", "stop_codon": "stop_codons"}.get(feature.lower())
            if destination:
                record[destination].append([a, b])
    if not records:
        raise ValueError("No GTF transcripts with gene_id/transcript_id were found")
    eligible, structural = defaultdict(list), defaultdict(list)
    transcript_exclusions = Counter()
    gene_reasons = defaultdict(set)
    gene_names = {}
    for record in records.values():
        gene_id = record["gene_id"]
        gene_names[gene_id] = record["gene_name"]
        try:
            model = segment_transcript(record, min_region_length)
        except ValueError as error:
            reason = str(error)
            transcript_exclusions[reason] += 1
            gene_reasons[gene_id].add(reason)
            continue
        structural[gene_id].append(model)
        reason = None
        if available_chrom_sizes is not None:
            if model["chrom"] not in available_chrom_sizes:
                reason = "contig_not_shared_by_all_tracks"
            elif max(end for start, end in model["exons"]) > available_chrom_sizes[model["chrom"]]:
                reason = "exons_outside_shared_contig_length"
        if reason:
            transcript_exclusions[reason] += 1
            gene_reasons[gene_id].add(reason)
        else:
            eligible[gene_id].append(model)

    def representative_order(model):
        return (-sum(model["lengths"].values()), -model["lengths"]["cds"],
                model["transcript_id"], model["chrom"], model["exons"], model["strand"])

    selected, selection_rows = [], []
    fallback_genes, genes_without_compatible_candidates = 0, 0
    for gene_id in sorted(gene_names):
        all_candidates = structural.get(gene_id, [])
        preferred = min(all_candidates, key=representative_order) if all_candidates else None
        candidates = eligible.get(gene_id, [])
        unavailable = len(all_candidates) - len(candidates)
        if candidates:
            model = min(candidates, key=representative_order)
            fallback = model is not preferred
            fallback_genes += fallback
            selected.append(model)
            selection_rows.append({
                "gene_id": gene_id, "gene_name": model["gene_name"], "status": "included",
                "reason": "available_contig_fallback" if fallback else "eligible_representative",
                "eligible_transcripts": len(candidates), "unavailable_candidates": unavailable,
                "preferred_transcript_id_before_contig_filter": preferred["transcript_id"],
                "preferred_chrom_before_contig_filter": preferred["chrom"],
                "selection_fallback": fallback,
                "transcript_id": model["transcript_id"], "chrom": model["chrom"],
                "strand": model["strand"], **model["lengths"],
            })
        else:
            genes_without_compatible_candidates += bool(all_candidates)
            selection_rows.append({
                "gene_id": gene_id, "gene_name": gene_names[gene_id], "status": "excluded",
                "reason": ";".join(sorted(gene_reasons[gene_id])), "eligible_transcripts": 0,
                "unavailable_candidates": unavailable,
                "preferred_transcript_id_before_contig_filter": preferred["transcript_id"] if preferred else "",
                "preferred_chrom_before_contig_filter": preferred["chrom"] if preferred else "",
                "selection_fallback": False,
                "transcript_id": "", "chrom": "", "strand": "", "utr5": 0, "cds": 0, "utr3": 0,
            })
    if not selected:
        raise ValueError("No compatible coding transcript has valid CDS and both UTRs; check annotation and BigWig contigs")
    selected.sort(key=lambda m: (m["chrom"], m["exons"][0][0], m["gene_id"]))
    summary = {
        "selection_policy": (("First restrict candidates to common BigWig contigs with consistent lengths and contained exons. "
                              if available_chrom_sizes is not None else "") + SELECTION_POLICY),
        "min_region_length": min_region_length,
        "coordinates": "0-based half-open genomic blocks, RNA 5prime-to-3prime order",
        "transcripts_seen": len(records), "genes_seen": len(gene_names),
        "transcripts_structurally_eligible": sum(map(len, structural.values())),
        "transcripts_eligible": sum(map(len, eligible.values())),
        "contig_filter_applied_before_selection": available_chrom_sizes is not None,
        "genes_with_available_contig_fallback": fallback_genes,
        "genes_no_compatible_contig_candidate": genes_without_compatible_candidates,
        "transcript_exclusions": dict(transcript_exclusions),
        "genes_included": len(selected), "genes_excluded": len(gene_names) - len(selected),
        "selected_with_explicit_stop_codon": sum(m["stop_codon_annotated"] for m in selected),
        "parse_counts": dict(parse_counts),
        "region_median_lengths": {region: statistics.median(m["lengths"][region] for m in selected)
                                  for region in REGIONS},
    }
    return selected, selection_rows, summary


def common_chromosome_sizes(headers):
    """Intersect all track contigs, rejecting disagreeing shared contig lengths."""
    if not headers:
        raise ValueError("At least one BigWig header is required for the common cohort")
    common = set(headers[0])
    union = set(headers[0])
    for header in headers[1:]:
        common.intersection_update(header)
        union.update(header)
    if not common:
        raise ValueError("Metaplot BigWigs share no common contigs")
    conflicts = [chrom for chrom in sorted(common)
                 if len({int(header[chrom]) for header in headers}) != 1]
    if conflicts:
        raise ValueError("Metaplot BigWigs disagree on contig lengths: " + ", ".join(conflicts[:10]))
    sizes = {chrom: int(headers[0][chrom]) for chrom in sorted(common)}
    if any(length <= 0 for length in sizes.values()):
        raise ValueError("BigWig contig lengths must be positive")
    return sizes, {"track_headers_checked": len(headers), "shared_contigs": len(common),
                   "union_contigs": len(union), "contigs_not_shared": sorted(union - common),
                   "shared_contig_lengths_consistent": True}


def scaled_region_profile(bigwig, chrom, blocks, strand, bins):
    """Exact interval-area means for fractional bins on concatenated RNA blocks.

    Missing BigWig coverage is zero. Memory is O(number of intervals in one
    transcript region + bins), never O(genome/transcript length).
    """
    import numpy as np

    intervals = []
    length = 0
    for start, end in blocks:
        for a, b, value in bigwig.intervals(chrom, start, end) or ():
            left, right = max(a, start), min(b, end)
            if left >= right or not np.isfinite(value):
                continue
            if strand == "+":
                intervals.append((length + left - start, length + right - start, value))
            else:
                intervals.append((length + end - right, length + end - left, value))
        length += end - start
    if length <= 0 or bins < 1:
        raise ValueError("Region length and bins must be positive")
    if not intervals:
        return np.zeros(bins, dtype=float)
    data = np.asarray(sorted(intervals), dtype=float)
    starts, ends, values = data.T
    areas = (ends - starts) * values
    prefix = np.concatenate(([0.0], np.cumsum(areas)))
    edges = np.linspace(0, length, bins + 1)
    index = np.searchsorted(ends, edges, side="right")
    integrals = prefix[index]
    valid = index < len(ends)
    j = index[valid]
    integrals[valid] += np.maximum(0, edges[valid] - starts[j]) * values[j]
    return np.diff(integrals) / (length / bins)


def load_models(path):
    with open(path) as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)
