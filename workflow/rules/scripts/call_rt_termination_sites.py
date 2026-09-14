"""Call sequence-unbiased RT peaks from the shared counted R2 signal.

BED6 outputs use RNA strand and 0-based crosslink coordinates. A physical
crosslink is not assigned to a motif or assumed to be the modified nucleotide.
R2-to-crosslink coordinates come from the common CIGAR-aware extractor.
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

if "snakemake" in globals():
    sys.path.insert(0, str(snakemake.params.script_dir))
from motif_utils import locus_key, revcomp, text_open


def load_counts(path):
    counts = defaultdict(dict)
    with text_open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, count, strand = line.rstrip().split("\t")[:6]
            pos, count = int(start), int(count)
            if int(end) != pos + 1 or pos < 0 or strand not in {"+", "-"} or count < 1:
                raise ValueError(f"Invalid counted RNA-strand 1-bp BED row: {line.strip()}")
            if pos in counts[(chrom, strand)]:
                raise ValueError(f"Duplicate R2 count key: {chrom}:{pos}:{strand}")
            counts[(chrom, strand)][pos] = count
    return counts


def select_rt_sites(counts, min_reads=40, neighbor_fold=1.5):
    """Published strict count > min_reads, and >= fold over either neighbor.

    Neighbors are immediately adjacent *genomic* positions on the same RNA
    strand, matching the reported pileup criterion. No sequence is inspected.
    """
    if min_reads < 0 or neighbor_fold <= 1:
        raise ValueError("min_reads >= 0 and neighbor_fold > 1 required")
    selected = {}
    for (chrom, strand), positions in counts.items():
        for pos, count in positions.items():
            left, right = positions.get(pos - 1, 0), positions.get(pos + 1, 0)
            if count > min_reads and count >= neighbor_fold * max(left, right):
                selected[(chrom, pos, strand)] = (count, left, right)
    return selected


def stream_selected_sites(path, min_reads=40, neighbor_fold=1.5):
    """Call same-strand adjacent-position peaks with O(selected sites) memory.

    Shared count BED is coordinate-sorted, with each chromosome contiguous.
    Only the preceding position on each strand is needed to assess neighbors.
    """
    if min_reads < 0 or neighbor_fold <= 1:
        raise ValueError("min_reads >= 0 and neighbor_fold > 1 required")
    selected, pending, seen_chroms = {}, {}, set()
    stats = Counter()
    current_chrom, last_position = None, -1

    def finalize(chrom, strand, candidate, right=0):
        pos, count, left = candidate
        if count > min_reads and count >= neighbor_fold * max(left, right):
            selected[(chrom, pos, strand)] = count, left, right

    with text_open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, count, strand = line.rstrip().split("\t")[:6]
            pos, count = int(start), int(count)
            if int(end) != pos + 1 or pos < 0 or count < 1 or strand not in {"+", "-"}:
                raise ValueError(f"Invalid counted BED row: {line.strip()}")
            if chrom != current_chrom:
                if chrom in seen_chroms:
                    raise ValueError("Count BED chromosome blocks must be contiguous")
                for previous_strand, candidate in pending.items():
                    finalize(current_chrom, previous_strand, candidate)
                current_chrom, last_position, pending = chrom, -1, {}
                seen_chroms.add(chrom)
            if pos < last_position:
                raise ValueError("Count BED must be coordinate sorted")
            previous = pending.get(strand)
            if previous is not None and pos == previous[0]:
                raise ValueError(f"Duplicate R2 count key: {chrom}:{pos}:{strand}")
            adjacent = previous is not None and previous[0] + 1 == pos
            if previous is not None:
                finalize(chrom, strand, previous, count if adjacent else 0)
            pending[strand] = pos, count, previous[1] if adjacent else 0
            last_position = pos
            stats["usable_read2"] += count
            stats["distinct_r2_positions"] += 1
    for strand, candidate in pending.items():
        finalize(current_chrom, strand, candidate)
    stats["selected_r2_positions"] = len(selected)
    return selected, stats


def mapped_crosslinks(mapping_path, selected, site_shift=1):
    """Preserve splice branches then aggregate unique physical crosslink sites."""
    sites = defaultdict(lambda: {"count": 0, "r2_positions": set(), "left": 0, "right": 0})
    mapping_totals = Counter()
    stats = Counter()
    with text_open(mapping_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "r2_pos0", "crosslink_pos0", "rna_strand", "count"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"R2 mapping must have header columns {sorted(required)}")
        for row in reader:
            chrom, strand = row["chrom"], row["rna_strand"]
            r2, crosslink, count = int(row["r2_pos0"]), int(row["crosslink_pos0"]), int(row["count"])
            key = chrom, r2, strand
            if key not in selected:
                continue
            if strand not in {"+", "-"} or crosslink < 0 or count < 1:
                raise ValueError(f"Invalid R2 mapping: {row}")
            mapping_totals[key] += count
            record = sites[(chrom, crosslink, strand)]
            record["count"] += count
            record["r2_positions"].add(r2)
            record["left"] = max(record["left"], selected[key][1])
            record["right"] = max(record["right"], selected[key][2])
            stats["selected_mapping_branches"] += 1
            stats["crosslink_supporting_reads"] += count
            if (r2 - crosslink) * (1 if strand == "+" else -1) != site_shift:
                stats["splice_adjusted_branches"] += 1
                stats["splice_adjusted_reads"] += count
    for key, (expected, _, _) in selected.items():
        observed = mapping_totals[key]
        if observed > expected:
            raise ValueError(f"Mapping support exceeds R2 count at {key}: {observed}>{expected}")
        stats["selected_r2_reads_without_valid_crosslink"] += expected - observed
    return sites, stats


def run(smk):
    import pysam
    selected, count_stats = stream_selected_sites(smk.input.counts, int(smk.params.min_reads), float(smk.params.neighbor_fold))
    sites, stats = mapped_crosslinks(smk.input.mapping, selected, int(smk.params.site_shift))
    stats.update(count_stats)
    sample, flank = str(smk.wildcards.sample), int(smk.params.flank)
    center_on = str(smk.params.sequence_center)
    if flank < 1 or center_on not in {"crosslink", "r2_start"}:
        raise ValueError("sequence_flank >= 1 and sequence_center crosslink/r2_start required")
    Path(str(smk.output.bed)).parent.mkdir(parents=True, exist_ok=True)
    stats["rt_sites"] = len(sites)
    with pysam.FastaFile(str(smk.input.fasta)) as fasta, open(str(smk.output.bed), "w") as bed, open(str(smk.output.fasta), "w") as out:
        lengths = dict(zip(fasta.references, fasta.lengths))
        rank = {chrom: i for i, chrom in enumerate(fasta.references)}
        for index, ((chrom, crosslink, strand), site) in enumerate(sorted(sites.items(), key=lambda item: (rank.get(item[0][0], len(rank)), item[0][1], item[0][2])), 1):
            name = f"{sample}_RT_{index}"
            starts = sorted(site["r2_positions"])
            bed.write(f"{chrom}\t{crosslink}\t{crosslink + 1}\t{name}\t{site['count']}\t{strand}\t{','.join(map(str, starts))}\t{site['left']}\t{site['right']}\n")
            # A unique physical site has one sequence; ambiguous R2-center sites
            # are omitted from discovery instead of duplicating sequence records.
            if center_on == "r2_start" and len(starts) != 1:
                stats["sequences_ambiguous_r2_center"] += 1
                continue
            center = crosslink if center_on == "crosslink" else starts[0]
            if chrom not in lengths or center - flank < 0 or center + flank + 1 > lengths[chrom]:
                stats["sequences_boundary_skipped"] += 1
                continue
            seq = fasta.fetch(chrom, center - flank, center + flank + 1).upper()
            if strand == "-":
                seq = revcomp(seq)
            if set(seq) - set("ACGT"):
                stats["sequences_ambiguous_bases_skipped"] += 1
                continue
            out.write(f">{name}|locus={locus_key(chrom, crosslink, strand)}|center={center_on}:{center}|count={site['count']}\n{seq}\n")
            stats["sequences"] += 1
    for key in ("usable_read2", "distinct_r2_positions", "selected_r2_positions", "rt_sites",
                "selected_mapping_branches", "crosslink_supporting_reads", "splice_adjusted_branches",
                "splice_adjusted_reads", "selected_r2_reads_without_valid_crosslink", "sequences",
                "sequences_ambiguous_r2_center", "sequences_boundary_skipped", "sequences_ambiguous_bases_skipped"):
        stats.setdefault(key, 0)
    report = {"sample": sample, "variant": str(smk.params.variant), **dict(stats), "min_reads_strictly_greater_than": int(smk.params.min_reads),
              "neighbor_fold_greater_or_equal": float(smk.params.neighbor_fold), "sequence_center": center_on,
              "sequence_flank": flank, "sequence_context": "genomic, RNA 5prime to 3prime; not spliced transcript windows",
              "selection": "R2 counts only; no motif anchoring", "unit": "distinct genomic crosslink locus and RNA strand"}
    with open(str(smk.output.stats), "w") as out:
        writer = csv.DictWriter(out, fieldnames=list(report), delimiter="\t")
        writer.writeheader()
        writer.writerow(report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
    with open(str(snakemake.log[0]), "w") as log, redirect_stdout(log), redirect_stderr(log):
        run(snakemake)
