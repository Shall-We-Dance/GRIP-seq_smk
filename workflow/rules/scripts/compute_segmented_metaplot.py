"""Stream transcript signal into segmented metaplot statistics and figures."""

from collections import Counter
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

sys.path.insert(0, str(snakemake.scriptdir))
from metagene import REGIONS, REGION_LABELS, load_models, scaled_region_profile


def moments_summary(total, total_square, count):
    if count == 0:
        return np.full_like(total, np.nan), np.full_like(total, np.nan)
    mean = total / count
    # Descriptive SEM across transcripts, not biological-replicate uncertainty.
    variance = np.maximum(0, (total_square - total * total / count) / max(1, count - 1))
    sem = np.sqrt(variance / count) if count > 1 else np.zeros_like(total)
    return mean, sem


def decorate_axis(ax, widths):
    edges = np.concatenate(([0], np.cumsum(widths)))
    for i, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        ax.axvspan(left, right, color=("#ddeaf4", "#ededed", "#f6e5d9")[i], alpha=0.35, zorder=-5)
        ax.text((left + right) / 2, 1.01, REGION_LABELS[i], transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9)
    for edge in edges[1:-1]:
        ax.axvline(edge, color="0.55", linestyle="--", linewidth=0.7)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xticks(edges)
    ax.set_xticklabels(["5′\nend", "CDS\nstart", "CDS\nend", "3′\nend"], fontsize=8)
    ax.grid(axis="y", alpha=0.18)


sample = str(snakemake.wildcards.sample)
bins = [int(snakemake.params.bins[region]) for region in REGIONS]
track_names = list(snakemake.params.tracks)
track_specs = list(snakemake.params.track_specs)
with open(snakemake.input.annotation_summary) as handle:
    annotation_summary = json.load(handle)
axis_scale = str(snakemake.params.axis_scale)
if axis_scale == "median_length":
    medians = annotation_summary["region_median_lengths"]
    widths = [medians[region] / medians["cds"] for region in REGIONS]
else:
    widths = [float(x) for x in snakemake.params.axis_widths]
if len(widths) != 3 or any(width <= 0 for width in widths):
    raise ValueError("Exactly three positive metaplot axis widths are required")
bin_widths = np.concatenate([np.full(n, width / n) for n, width in zip(bins, widths)])
x = np.cumsum(bin_widths) - bin_widths / 2
n_bins = sum(bins)
accumulators = {
    track: {"sum": np.zeros(n_bins), "sum_square": np.zeros(n_bins),
            "shape_sum": np.zeros(n_bins), "shape_sum_square": np.zeros(n_bins),
            "mass": np.zeros(n_bins), "nonzero": 0}
    for track in track_names
}
counts = Counter(annotation_genes=0, included_genes=0, excluded_missing_contig=0, excluded_out_of_bounds=0)
for path in snakemake.output:
    Path(path).parent.mkdir(parents=True, exist_ok=True)

with ExitStack() as stack:
    gene_handle = stack.enter_context(open(snakemake.output.gene_signal, "w"))
    gene_writer = csv.DictWriter(gene_handle, fieldnames=["sample", "track", "gene_id", "transcript_id",
        "chrom", "strand", "mature_length", "utr5_length", "cds_length", "utr3_length",
        "signal_mass", "mean_cpm"], delimiter="\t")
    gene_writer.writeheader()
    paths = list(snakemake.input.bws)
    handles = [stack.enter_context(pyBigWig.open(str(path))) for path in paths]
    chrom_sizes = [handle.chroms() for handle in handles]
    for model in load_models(snakemake.input.models):
        counts["annotation_genes"] += 1
        chrom = model["chrom"]
        if any(chrom not in sizes for sizes in chrom_sizes):
            raise ValueError(f"Shared metaplot annotation contains unavailable contig {chrom}; rebuild the common cohort")
        end = max(b for a, b in model["exons"])
        if any(end > sizes[chrom] for sizes in chrom_sizes):
            raise ValueError(f"Shared metaplot annotation exceeds contig {chrom} bounds; rebuild the common cohort")
        counts["included_genes"] += 1
        lengths = [model["lengths"][region] for region in REGIONS]
        nt_widths = np.concatenate([np.full(n, length / n) for n, length in zip(bins, lengths)])
        for name, spec in zip(track_names, track_specs):
            index = spec["plus"] if model["strand"] == "+" else spec["minus"]
            bigwig = handles[index]
            profile = np.concatenate([
                scaled_region_profile(bigwig, chrom, model["regions"][region], model["strand"], n)
                for region, n in zip(REGIONS, bins)
            ])
            if np.any(profile < -1e-8):
                raise ValueError(f"{name}: metaplot density requires nonnegative coverage tracks")
            profile = np.maximum(profile, 0)
            mass = profile * nt_widths
            acc = accumulators[name]
            acc["sum"] += profile
            acc["sum_square"] += profile * profile
            acc["mass"] += mass
            total_mass = mass.sum()
            gene_writer.writerow({"sample": sample, "track": name, "gene_id": model["gene_id"],
                "transcript_id": model["transcript_id"], "chrom": chrom, "strand": model["strand"],
                "mature_length": sum(lengths), "utr5_length": lengths[0], "cds_length": lengths[1],
                "utr3_length": lengths[2], "signal_mass": float(total_mass),
                "mean_cpm": float(total_mass / sum(lengths))})
            if total_mass > 0:
                shape = profile / (total_mass / sum(lengths))
                acc["shape_sum"] += shape
                acc["shape_sum_square"] += shape * shape
                acc["nonzero"] += 1

if not counts["included_genes"]:
    raise ValueError("No annotated coding transcripts occur in every requested BigWig; check contig names")
for path in snakemake.output:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
statistics = {
    "sample": sample, "axis_scale": axis_scale, "axis_widths": dict(zip(REGIONS, widths)),
    "bins": dict(zip(REGIONS, bins)), "cohort": dict(counts),
    "missing_bigwig_signal": "zero; transcript candidates are restricted to contigs shared across all configured samples/tracks before representative selection",
    "selection_policy": annotation_summary["selection_policy"],
    "density_definition": "Sum exact signal area per scaled bin over selected transcripts; "
        "divide by total mapped transcript signal and display-bin width (area integrates to 1). "
        "Overlapping genes may count a genomic observation more than once. R2/crosslink tracks "
        "are read weighted; BAM coverage is covered-base weighted. No CLIPper cluster filtering.",
    "relative_density_definition": "Per-transcript CPM divided by its mean CPM over mature RNA; "
        "equal-weight mean over nonzero transcripts. Other CPM means include zero-signal genes.",
    "uncertainty": "SEM columns describe across-transcript variability, not biological replication.",
    "tracks": {},
}
fig, axes = plt.subplots(len(track_names), 3, figsize=(16, 3.25 * len(track_names)), squeeze=False)
columns = ["sample", "track", "region", "bin", "global_bin", "region_fraction", "x", "x_width",
           "probability_density", "mean_cpm", "sem_cpm", "n_transcripts", "relative_density",
           "relative_density_sem", "n_nonzero_transcripts"]
with open(snakemake.output.data, "w") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    for row_index, name in enumerate(track_names):
        acc = accumulators[name]
        mean, sem = moments_summary(acc["sum"], acc["sum_square"], counts["included_genes"])
        shape, shape_sem = moments_summary(acc["shape_sum"], acc["shape_sum_square"], acc["nonzero"])
        total_mass = acc["mass"].sum()
        density = (acc["mass"] / total_mass / bin_widths if total_mass > 0
                   else np.full(n_bins, np.nan))
        statistics["tracks"][name] = {
            "nonzero_genes": acc["nonzero"], "mapped_transcript_signal_mass": float(total_mass),
            "density_integral": float(np.sum(density * bin_widths)) if total_mass > 0 else None,
            "strand_aware": track_specs[row_index]["plus"] != track_specs[row_index]["minus"],
            "status": "ok" if total_mass > 0 else "no_signal",
        }
        index = 0
        for region, n in zip(REGIONS, bins):
            for bin_index in range(n):
                writer.writerow(dict(zip(columns, [sample, name, region, bin_index + 1, index + 1,
                    (bin_index + 0.5) / n, x[index], bin_widths[index], density[index],
                    mean[index], sem[index], counts["included_genes"], shape[index],
                    shape_sem[index], acc["nonzero"]])))
                index += 1
        for ax, data, label in zip(axes[row_index], (density, mean, shape),
                                  ("Signal probability density", "Mean CPM", "Equal-gene relative density")):
            ax.plot(x, data, color="#296997", linewidth=1.3)
            decorate_axis(ax, widths)
            ax.set_ylabel(label)
            ax.set_title(name + (" (RNA strand)" if statistics["tracks"][name]["strand_aware"] else " (unstranded)"),
                         fontsize=9, pad=25)
            if total_mass == 0:
                ax.text(0.5, 0.5, "No signal in selected transcripts", transform=ax.transAxes, ha="center")
fig.suptitle(f"{sample}: spliced 5′UTR–CDS–3′UTR; {counts['included_genes']:,} genes", fontsize=12)
fig.supxlabel("Each region normalized independently; displayed widths = " +
              ("median region length / median CDS length" if axis_scale == "median_length" else "configured ratios"), fontsize=10)
fig.tight_layout(rect=(0, 0.025, 1, 0.965))
fig.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close(fig)
with open(snakemake.output.stats, "w") as handle:
    json.dump(statistics, handle, indent=2, sort_keys=True)
    handle.write("\n")
