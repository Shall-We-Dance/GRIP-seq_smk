"""Compare sample-normalized densities by experimental IP/Input group.

Shading represents biological-sample SD when n>=2; a single Input has no
uncertainty band. The exported table also retains mean CPM and equal-gene
relative densities for each condition.
"""

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(snakemake.scriptdir))
from metagene import REGIONS, REGION_LABELS

METRICS = ("probability_density", "mean_cpm", "relative_density")
profiles = {}
coordinates = None
widths = None
for path in snakemake.input.profiles:
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["sample"], row["track"])
            profiles.setdefault(key, []).append(row)
for key, rows in profiles.items():
    current = np.array([float(row["x"]) for row in rows])
    current_widths = np.array([float(row["x_width"]) for row in rows])
    if coordinates is None:
        coordinates, widths = current, current_widths
        template_rows = rows
    elif (len(current) != len(coordinates) or not np.allclose(current, coordinates)
          or not np.allclose(current_widths, widths)):
        raise ValueError("Cannot combine profiles with different region normalization or bin counts")

sample_names = list(snakemake.params.samples)
track_names = list(snakemake.params.tracks)
configured_groups = dict(snakemake.params.groups)
groups = {
    name: {role: [sample for sample in names if sample in sample_names]
           for role, names in roles.items() if role in ("IP", "Input")}
    for name, roles in configured_groups.items()
}
groups = {name: roles for name, roles in groups.items() if any(roles.values())}
if not groups:
    # Without experimental metadata each sample gets its own curve, with no
    # invented replicate grouping or pooled all-sample mean.
    groups = {"Samples": {sample: [sample] for sample in sample_names}}
else:
    assigned = {sample for roles in groups.values() for names in roles.values() for sample in names}
    if set(sample_names) - assigned:
        groups["Ungrouped samples"] = {sample: [sample] for sample in sample_names if sample not in assigned}
region_widths = [sum(float(row["x_width"]) for row in template_rows if row["region"] == region)
                 for region in REGIONS]
edges = np.concatenate(([0], np.cumsum(region_widths)))
fig, axes = plt.subplots(len(track_names), len(groups), figsize=(5.0 * len(groups), 3.8 * len(track_names)),
                         squeeze=False, sharex=True)
for path in snakemake.output:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
fields = ["group", "role", "track", "metric", "region", "bin", "global_bin", "x", "x_width",
          "mean", "sd_samples", "n_samples"]
with open(snakemake.output.data, "w") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for column_index, (group, roles) in enumerate(groups.items()):
        for row_index, track in enumerate(track_names):
            ax = axes[row_index, column_index]
            for role_index, (role, names) in enumerate(roles.items()):
                if not names:
                    continue
                color = {"IP": "#bd473d", "Input": "#336e9f"}.get(role, plt.cm.tab20(role_index % 20))
                for metric in METRICS:
                    matrix = np.array([[float(row[metric]) for row in profiles[(name, track)]] for name in names])
                    valid_n = np.isfinite(matrix).sum(axis=0)
                    mean = np.divide(np.nansum(matrix, axis=0), valid_n,
                                     out=np.full(matrix.shape[1], np.nan), where=valid_n > 0)
                    differences = np.where(np.isfinite(matrix), matrix - mean, 0)
                    sd = np.sqrt(np.divide((differences * differences).sum(axis=0), valid_n - 1,
                                           out=np.full(matrix.shape[1], np.nan), where=valid_n > 1))
                    for i, row in enumerate(template_rows):
                        writer.writerow(dict(zip(fields, [group, role, track, metric, row["region"],
                            row["bin"], row["global_bin"], coordinates[i], widths[i], mean[i], sd[i], valid_n[i]])))
                    if metric == "probability_density":
                        for values in matrix:
                            ax.plot(coordinates, values, color=color, alpha=0.22, linewidth=0.6)
                        n_valid = np.isfinite(matrix).any(axis=1).sum()
                        ax.plot(coordinates, mean, color=color, linewidth=1.8, label=f"{role} (n={n_valid})")
                        if np.any(valid_n > 1):
                            ax.fill_between(coordinates, np.maximum(0, mean - sd), mean + sd,
                                            color=color, alpha=0.17, linewidth=0)
            for i, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
                ax.axvspan(left, right, color=("#ddeaf4", "#ededed", "#f6e5d9")[i], alpha=0.3, zorder=-5)
                ax.text((left + right) / 2, 1.01, REGION_LABELS[i], transform=ax.get_xaxis_transform(),
                        ha="center", fontsize=9)
            for edge in edges[1:-1]:
                ax.axvline(edge, color="0.55", linestyle="--", linewidth=0.7)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_xticks(edges)
            ax.set_xticklabels(["5′\nend", "CDS\nstart", "CDS\nend", "3′\nend"], fontsize=8)
            ax.set_title(f"{group}\n{track}", fontsize=9, pad=22)
            ax.set_ylabel("Signal probability density")
            ax.grid(axis="y", alpha=0.18)
            ax.legend(fontsize=8, frameon=False)
fig.suptitle("Spliced 5′UTR–CDS–3′UTR metaplots\nSample means; shading = sample SD (n ≥ 2)", fontsize=12)
fig.supxlabel("Independent region normalization; " + ("median-length segment widths" if
              str(snakemake.params.axis_scale) == "median_length" else "configured segment widths"), fontsize=10)
fig.tight_layout(rect=(0, 0.025, 1, 0.95))
fig.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close(fig)
