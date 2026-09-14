"""Sequence-unbiased site-level motif offsets and local shifted background enrichment.

One foreground site is paired to one deterministic nearby genomic position on
its RNA strand. Background centers are sampled without inspecting motif status.
Foreground/background discordance gives a one-sided exact paired-binomial test
(McNemar). BH spans every tested motif and offset within each sample. These are
exploratory locus-level tests, not tests over independent biological replicates;
neighboring genomic loci can be correlated. Read weights are descriptive only.
"""
import csv
from bisect import bisect_left
import json
import math
import sys
from collections import Counter, defaultdict
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

if "snakemake" in globals():
    sys.path.insert(0, str(snakemake.params.script_dir))
from motif_utils import (bh_adjust, is_holdout, locus_key, locus_partition, motif_matches, motif_regex,
                         parse_meme_consensuses, revcomp, stable_hash, text_open)


def read_sites(path, with_r2=False):
    sites = {}
    with text_open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, position, end, name, count, strand = fields[:6]
            position, count = int(position), int(count)
            if int(end) != position + 1 or position < 0 or count < 1 or strand not in {"+", "-"}:
                raise ValueError(f"Invalid crosslink BED row: {line.strip()}")
            key = (chrom, position, strand)
            if key in sites:
                raise ValueError(f"Crosslink BED must be unique by locus and RNA strand: {key}")
            if with_r2:
                starts = [int(value) for value in fields[6].split(",")] if len(fields) > 6 and fields[6] else []
                sites[key] = {"count": count, "r2_positions": starts}
            else:
                sites[key] = count
    return sites


def load_blacklist(paths):
    raw = defaultdict(list)
    for path in paths:
        with text_open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                chrom, start, end = line.split()[:3]
                start, end = int(start), int(end)
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid blacklist interval: {line.strip()}")
                raw[chrom].append((start, end))
    result = {}
    for chrom, intervals in raw.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = merged[-1][0], max(merged[-1][1], end)
            else:
                merged.append((start, end))
        result[chrom] = ([start for start, end in merged], [end for start, end in merged])
    return result


def overlaps_blacklist(blacklist, chrom, start, end):
    starts, ends = (blacklist or {}).get(chrom, ([], []))
    last_before_end = bisect_left(starts, end) - 1
    return last_before_end >= 0 and ends[last_before_end] > start


def sequence_at(fasta, chrom, position, strand, flank, chrom_lengths, blacklist=None):
    if chrom not in chrom_lengths or position - flank < 0 or position + flank + 1 > chrom_lengths[chrom]:
        return None, "boundary"
    if overlaps_blacklist(blacklist, chrom, position - flank, position + flank + 1):
        return None, "blacklist"
    sequence = fasta.fetch(chrom, position - flank, position + flank + 1).upper()
    if set(sequence) - set("ACGT"):
        return None, "ambiguous_bases"
    return (revcomp(sequence) if strand == "-" else sequence), None


def matched_background(fasta, chrom, position, strand, flank, lengths, min_distance, max_distance, seed, max_attempts, blacklist=None, diagnostics=None, block_size=0):
    key = locus_key(chrom, position, strand)
    for attempt in range(max_attempts):
        draw = stable_hash(f"background|{key}|{attempt}", seed)
        sign = 1 if draw & 1 else -1
        distance = min_distance + (draw >> 1) % (max_distance - min_distance + 1)
        background_pos = position + sign * distance
        if block_size and (background_pos - flank < (position // block_size) * block_size or background_pos + flank >= (position // block_size + 1) * block_size):
            if diagnostics is not None:
                diagnostics["background_candidates_cross_holdout_block_rejected"] += 1
            continue
        sequence, reason = sequence_at(fasta, chrom, background_pos, strand, flank, lengths, blacklist)
        if diagnostics is not None:
            diagnostics["background_candidates_examined"] += 1
            if reason is not None:
                diagnostics["background_candidates_" + reason + "_rejected"] += 1
        if sequence is not None:
            return sequence, background_pos
    return None, None


def new_accumulator(window):
    size = 2 * window + 1
    return {"sites": 0, "reads": 0, "foreground_any": 0, "background_any": 0,
            "foreground_multiple": 0, "background_multiple": 0,
            "foreground": [0] * size, "background": [0] * size,
            "foreground_reads": [0] * size, "background_reads": [0] * size,
            "foreground_only": [0] * size, "background_only": [0] * size}


def add_pair(acc, foreground_offsets, background_offsets, read_count, window):
    acc["sites"] += 1
    acc["reads"] += read_count
    acc["foreground_any"] += bool(foreground_offsets)
    acc["background_any"] += bool(background_offsets)
    acc["foreground_multiple"] += len(foreground_offsets) > 1
    acc["background_multiple"] += len(background_offsets) > 1
    for offsets, count_key, read_key in ((foreground_offsets, "foreground", "foreground_reads"),
                                         (background_offsets, "background", "background_reads")):
        for offset in offsets:
            index = offset + window
            acc[count_key][index] += 1
            acc[read_key][index] += read_count
    for offset in foreground_offsets - background_offsets:
        acc["foreground_only"][offset + window] += 1
    for offset in background_offsets - foreground_offsets:
        acc["background_only"][offset + window] += 1


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def make_rows(name, motif, acc, window, site_shift=1, r2_subset=False):
    from scipy.stats import binomtest
    rows = []
    n_sites, n_reads = acc["sites"], acc["reads"]
    for index, offset in enumerate(range(-window, window + 1)):
        foreground, background = acc["foreground"][index], acc["background"][index]
        fg_only, bg_only = acc["foreground_only"][index], acc["background_only"][index]
        discordant = fg_only + bg_only
        pvalue = (binomtest(fg_only, discordant, 0.5, alternative="greater").pvalue if discordant else 1.0) if n_sites else None
        rows.append({"motif": name, "iupac": motif["sequence"], "anchor_pos_1based": motif["anchor_pos"],
                     "source": motif["source"], "motif_anchor_offset_from_crosslink_rna": offset,
                     "crosslink_offset_from_motif_anchor_rna": -offset,
                     "r2_start_offset_from_motif_anchor_rna": site_shift - offset if r2_subset else None,
                     "coordinate_population": "unique_contiguous_R2" if r2_subset else "all_eligible_crosslinks",
                     "effective_site_pairs": n_sites, "supporting_reads": n_reads,
                     "foreground_sites": foreground, "background_sites": background,
                     "foreground_site_fraction": ratio(foreground, n_sites),
                     "background_site_fraction": ratio(background, n_sites),
                     "foreground_read_weight": acc["foreground_reads"][index],
                     "background_read_weight": acc["background_reads"][index],
                     "foreground_read_fraction": ratio(acc["foreground_reads"][index], n_reads),
                     "background_read_fraction": ratio(acc["background_reads"][index], n_reads),
                     "enrichment_ratio_pseudocount_0p5": (foreground + 0.5) / (background + 0.5) if n_sites else None,
                     "foreground_only_pairs": fg_only, "background_only_pairs": bg_only,
                     "discordant_pairs": discordant, "pvalue_one_sided_paired": pvalue})
    return rows


def summarize_motif(motif, acc, rows):
    n = acc["sites"]
    if n and max(acc["foreground"], default=0) > 0:
        max_count = max(acc["foreground"])
        peak_offsets = [row["motif_anchor_offset_from_crosslink_rna"] for row in rows if row["foreground_sites"] == max_count]
        max_excess = max(row["foreground_sites"] - row["background_sites"] for row in rows)
        excess_offsets = [row["motif_anchor_offset_from_crosslink_rna"] for row in rows if row["foreground_sites"] - row["background_sites"] == max_excess]
        significant = [row["motif_anchor_offset_from_crosslink_rna"] for row in rows if row["qvalue_bh_sample"] is not None and row["qvalue_bh_sample"] < 0.05 and row["foreground_sites"] > row["background_sites"]]
    else:
        peak_offsets, excess_offsets, significant = [], [], []
    return {**motif, "effective_site_pairs": n, "supporting_reads": acc["reads"],
            "sites_with_any_foreground_match": acc["foreground_any"],
            "sites_with_any_background_match": acc["background_any"],
            "sites_with_multiple_foreground_offsets": acc["foreground_multiple"],
            "sites_with_multiple_background_offsets": acc["background_multiple"],
            "peak_frequency_offsets_nt": peak_offsets, "peak_excess_over_background_offsets_nt": excess_offsets,
            "significant_enrichment_offsets_q_lt_0p05": significant,
            "offset_estimate_is_unique": len(peak_offsets) == 1,
            "physical_crosslink_coordinates_modified": False,
            "status": "no_effective_sites" if n == 0 else "no_motif_matches" if not peak_offsets else "completed"}


def plot_distribution(rows, summaries, r2_rows, r2_summaries, sample, out_png, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    fig, axes = plt.subplots(max(1, len(summaries)), 2, squeeze=False,
                             figsize=(12, max(3.6, 2.8 * len(summaries))))
    for index, summary in enumerate(summaries):
        for col, (xkey, xlabel) in enumerate((("motif_anchor_offset_from_crosslink_rna", "Motif anchor relative to crosslink (nt; RNA 3′ positive)"),
                                             ("r2_start_offset_from_motif_anchor_rna", "R2 start relative to motif anchor (nt; unique contiguous R2 subset)"))):
            ax = axes[index, col]
            panel_rows = rows if col == 0 else r2_rows
            panel_summary = summary if col == 0 else r2_summaries[index]
            motif_rows = [row for row in panel_rows if row["motif"] == summary["name"]]
            ordered = sorted(motif_rows, key=lambda row: row[xkey])
            x = [row[xkey] for row in ordered]
            if panel_summary["effective_site_pairs"]:
                ax.plot(x, [row["foreground_site_fraction"] for row in ordered], label="Distinct sites", color="#1863ac", lw=1.6)
                ax.plot(x, [row["background_site_fraction"] for row in ordered], label="Local shifted background", color="#777777", lw=1.2)
                ax.plot(x, [row["foreground_read_fraction"] for row in ordered], label="Read weighted (descriptive)", color="#e38c32", lw=1, alpha=.75)
                significant = [row for row in ordered if row["qvalue_bh_sample"] is not None and row["qvalue_bh_sample"] < .05 and row["foreground_sites"] > row["background_sites"]]
                if significant:
                    ax.scatter([row[xkey] for row in significant], [row["foreground_site_fraction"] for row in significant], s=13, color="#1863ac", marker="*", zorder=4)
            else:
                ax.text(.5, .5, "No eligible foreground/background pairs", ha="center", va="center", transform=ax.transAxes)
            ax.axvline(0, color="#bbbbbb", lw=.8, ls="--")
            ax.set(xlabel=xlabel, ylabel="Fraction", title=f"{summary['name']} ({summary['sequence']}; anchor {summary['anchor_pos']}) · n={panel_summary['effective_site_pairs']:,}")
            ax.set_ylim(bottom=0)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.spines[["top", "right"]].set_visible(False)
            if index == 0 and col == 0 and summary["effective_site_pairs"]:
                ax.legend(frameon=False, fontsize=8)
    fig.suptitle(sample + " · RNA-oriented motif distributions", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.45 / fig.get_figheight()))
    fig.savefig(out_png, dpi=180)
    fig.savefig(out_pdf)
    plt.close(fig)


def run(smk):
    import pysam
    sample = str(smk.wildcards.sample)
    window, seed = int(smk.params.window), int(smk.params.seed)
    holdout_fraction, site_shift = float(smk.params.holdout_fraction), int(smk.params.site_shift)
    block_size, holdout_guard = int(smk.params.holdout_block_size), int(smk.params.holdout_guard)
    motifs = [dict(m) for m in smk.params.motifs]
    for motif in motifs:
        motif.setdefault("source", "predefined")
        motif_regex(motif["sequence"])
        if not 1 <= int(motif["anchor_pos"]) <= len(motif["sequence"]):
            raise ValueError(f"Invalid anchor position for {motif}")
    if list(smk.input.meme):
        motifs.extend(parse_meme_consensuses(str(smk.input.meme[0]), float(smk.params.consensus_mass)))
    if len(set(m["name"] for m in motifs)) != len(motifs):
        raise ValueError("Motif names must be unique")
    if not motifs or window < 0 or not 0 <= holdout_fraction < 1:
        raise ValueError("At least one motif, window >=0, and holdout_fraction in [0,1) required")
    flank = window + max(len(m["sequence"]) - 1 for m in motifs)
    if holdout_guard < flank or block_size <= 2 * holdout_guard:
        raise ValueError("Holdout guard must cover every scan window, with block_size > 2*guard")
    minimum, maximum = int(smk.params.background_min_distance), int(smk.params.background_max_distance)
    if minimum <= 2 * flank or maximum < minimum or int(smk.params.background_attempts) < 1:
        raise ValueError("Background min distance must exceed twice the scan flank; max >= min and attempts >=1")
    sites = read_sites(smk.input.bed, with_r2=True)
    blacklist = load_blacklist(getattr(smk.input, "blacklist", []))
    stats = Counter({key: 0 for key in ("matched_site_pairs", "matched_supporting_reads", "matched_holdout_sites",
        "matched_holdout_guard_sites", "foreground_boundary_skipped", "foreground_ambiguous_bases_skipped",
        "foreground_blacklist_skipped", "no_valid_local_background_skipped", "background_candidates_examined",
        "background_candidates_boundary_rejected", "background_candidates_ambiguous_bases_rejected",
        "background_candidates_blacklist_rejected", "background_candidates_cross_holdout_block_rejected",
        "r2_panel_eligible_site_pairs", "r2_panel_spliced_ambiguous_or_missing_mapping_excluded")})
    stats.update(input_sites=len(sites), input_supporting_reads=sum(site["count"] for site in sites.values()))
    accumulators = {motif["name"]: new_accumulator(window) for motif in motifs}
    r2_accumulators = {motif["name"]: new_accumulator(window) for motif in motifs}
    Path(str(smk.output.table)).parent.mkdir(parents=True, exist_ok=True)
    with pysam.FastaFile(str(smk.input.fasta)) as fasta, open(str(smk.output.pairs), "w") as pairs:
        lengths = dict(zip(fasta.references, fasta.lengths))
        pairs.write("chrom\tcrosslink_pos0\trna_strand\tsupporting_reads\tbackground_pos0\tgenomic_background_distance\tdiscovery_holdout\tholdout_partition\tunique_contiguous_r2\n")
        for (chrom, pos, strand), site in sites.items():
            count = site["count"]
            key = locus_key(chrom, pos, strand)
            partition = locus_partition(key, holdout_fraction, seed, block_size, holdout_guard)
            r2_starts = site["r2_positions"]
            r2_eligible = len(r2_starts) == 1 and (r2_starts[0] - pos) * (1 if strand == "+" else -1) == site_shift
            if not r2_eligible:
                stats["r2_panel_spliced_ambiguous_or_missing_mapping_excluded"] += 1
            foreground, reason = sequence_at(fasta, chrom, pos, strand, flank, lengths, blacklist)
            if foreground is None:
                stats["foreground_" + reason + "_skipped"] += 1
                continue
            background, bg_pos = matched_background(fasta, chrom, pos, strand, flank, lengths, minimum, maximum,
                                                   seed, int(smk.params.background_attempts), blacklist, stats, block_size)
            if background is None:
                stats["no_valid_local_background_skipped"] += 1
                continue
            heldout = partition == "holdout"
            stats["matched_site_pairs"] += 1
            stats["matched_supporting_reads"] += count
            stats["matched_holdout_sites"] += heldout
            stats["matched_holdout_guard_sites"] += partition == "boundary_guard"
            stats["r2_panel_eligible_site_pairs"] += r2_eligible
            pairs.write(f"{chrom}\t{pos}\t{strand}\t{count}\t{bg_pos}\t{bg_pos - pos}\t{int(heldout)}\t{partition}\t{int(r2_eligible)}\n")
            for motif in motifs:
                if motif["source"] == "de_novo_consensus_holdout" and not heldout:
                    continue
                fg_offsets = motif_matches(foreground, flank, motif["sequence"], int(motif["anchor_pos"]), window)
                bg_offsets = motif_matches(background, flank, motif["sequence"], int(motif["anchor_pos"]), window)
                add_pair(accumulators[motif["name"]], fg_offsets, bg_offsets, count, window)
                if r2_eligible:
                    add_pair(r2_accumulators[motif["name"]], fg_offsets, bg_offsets, count, window)
    rows, r2_rows = [], []
    for motif in motifs:
        rows.extend(make_rows(motif["name"], motif, accumulators[motif["name"]], window, site_shift))
        r2_rows.extend(make_rows(motif["name"], motif, r2_accumulators[motif["name"]], window, site_shift, r2_subset=True))
    all_rows = rows + r2_rows
    qvalues = bh_adjust([row["pvalue_one_sided_paired"] for row in all_rows])
    for row, qvalue in zip(all_rows, qvalues):
        row["qvalue_bh_sample"] = qvalue
    summaries = [summarize_motif(m, accumulators[m["name"]], [row for row in rows if row["motif"] == m["name"]]) for m in motifs]
    r2_summaries = [summarize_motif(m, r2_accumulators[m["name"]], [row for row in r2_rows if row["motif"] == m["name"]]) for m in motifs]
    for filename, table_rows in ((smk.output.table, rows), (smk.output.r2_table, r2_rows)):
        with open(str(filename), "w") as out:
            writer = csv.DictWriter(out, fieldnames=list(table_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows({key: "NA" if value is None else value for key, value in row.items()} for row in table_rows)
    report = {"sample": sample, "variant": str(getattr(smk.params, "variant", "noblacklist")), "counts": dict(stats), "motifs": summaries, "r2_contiguous_subset_motifs": r2_summaries, "window_nt": window,
              "site_shift_upstream_rna_nt": site_shift, "seed": seed, "de_novo_holdout_fraction": holdout_fraction,
              "de_novo_holdout_block_size": block_size, "de_novo_holdout_guard": holdout_guard,
              "r2_coordinate_population": "only unique R2 centers whose true genomic RNA-oriented distance equals site_shift; splice/ambiguous mappings are excluded",
              "background": {"method": "one deterministic nearby genomic center per foreground site, same chromosome and RNA strand",
                             "min_distance_nt": minimum, "max_distance_nt": maximum,
                             "motif_status_used_to_select_background": False,
                             "foreground_and_background_windows_exclude_blacklist": bool(blacklist),
                             "background_overlaps_own_foreground_window": False,
                             "background_whole_window_in_foreground_genomic_holdout_block": True},
              "offset_peak_inference": "exploratory; data-selected maxima, not pre-registered single-offset tests",
              "statistical_test": "one-sided exact paired binomial (McNemar) on discordant foreground/background site pairs",
              "multiple_testing": "Benjamini-Hochberg across all tested motif/offset combinations and both crosslink/R2 populations within this sample",
              "zero_site_behavior": "zero counts; fractions, enrichment, p and q are NA; no offset estimate",
              "read_weights": "descriptive; supporting reads are not independent site or biological replicates",
              "sequence_context": "genomic windows in RNA orientation; no transcript splicing or exon matching",
              "caveats": ["Nearby loci and overlapping windows may be correlated; locus-level p/q values are exploratory and do not replace biological-replicate inference.",
                          "Local genomic backgrounds control location and strand, but not transcript isoform, exon membership, expression, or sequence composition exactly.",
                          "All overlapping motif occurrences are counted once per site at each offset; per-offset fractions need not sum to one.",
                          "De novo consensuses use a declared PWM probability-mass IUPAC approximation, not FIMO PWM p-values, and use held-out loci only.",
                          "For even-width de novo motifs the left central base is the anchor; this is a plotting convention, not a modified-nucleotide assignment.",
                          "R2-relative coordinates use verified unique contiguous R2 mappings only; splice/multiple/missing R2 mappings have a separate exclusion count and remain in the crosslink population.",
                          "Peak estimates preserve every tied maximum and never change physical crosslink coordinates."]}
    Path(str(smk.output.stats)).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    plot_distribution(rows, summaries, r2_rows, r2_summaries, sample, str(smk.output.png), str(smk.output.pdf))
    print(json.dumps({"sample": sample, "counts": dict(stats), "motifs": summaries}, indent=2))


if __name__ == "__main__":
    Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
    with open(str(snakemake.log[0]), "w") as log, redirect_stdout(log), redirect_stderr(log):
        run(snakemake)
