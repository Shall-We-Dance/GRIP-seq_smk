"""Combine small deepTools profile tables into a readable two-panel figure."""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


sample_names = list(snakemake.params.samples)
n_bins = int(snakemake.params.n_bins)
upstream = int(snakemake.params.upstream)
body = int(snakemake.params.body)
downstream = int(snakemake.params.downstream)
bin_size = int(snakemake.params.bin_size)


def profiles_from_file(path):
    profiles = []
    with open(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            # deepTools writes two header rows ("bin labels" and "bins")
            # followed by: sample, group, value1 ... valueN.
            if not row or row[0] in ("bin labels", "bins") or row[0].startswith("#"):
                continue
            numeric = []
            for value in row[2:]:
                try:
                    numeric.append(float(value))
                except ValueError:
                    pass
            if len(numeric) == n_bins:
                profiles.append(np.asarray(numeric, dtype=float))
    if len(profiles) < 2:
        raise ValueError(f"Expected two profiles in {path}, found {len(profiles)}")
    return profiles[:2]


fragment_profiles = []
termination_profiles = []
for path in snakemake.input:
    profiles = profiles_from_file(str(path))
    fragment_profiles.append(profiles[0])
    termination_profiles.append(profiles[1])

x = np.arange(n_bins)
tss = upstream // bin_size
tes = (upstream + body) // bin_size
fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
titles = ["Mapped-fragment coverage (BAM CPM BigWig)", "Read 2 first-base signal (1-bp CPM BigWig)"]
for ax, data, title in zip(axes, (fragment_profiles, termination_profiles), titles):
    matrix = np.vstack(data)
    colors = plt.cm.tab20(np.linspace(0, 1, max(20, len(sample_names))))
    for index, (sample, profile) in enumerate(zip(sample_names, matrix)):
        ax.plot(x, profile, lw=0.8, alpha=0.55, color=colors[index], label=sample)
    ax.plot(x, np.nanmean(matrix, axis=0), color="black", lw=2.2, label="all-sample mean")
    ax.axvline(tss, color="0.45", ls="--", lw=0.8)
    ax.axvline(tes, color="0.45", ls="--", lw=0.8)
    ax.set_title(title)
    ax.set_ylabel("Mean CPM")
    ax.grid(axis="y", alpha=0.2)

axes[-1].set_xticks([0, tss, tes, n_bins - 1])
axes[-1].set_xticklabels([f"-{upstream / 1000:g} kb", "TSS", "TES", f"+{downstream / 1000:g} kb"])
axes[-1].set_xlabel("Scaled gene body (RNA 5′ → 3′)")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.815, 0.5), fontsize=7)
fig.suptitle("GRIP-seq TSS–TES metaplots: all experiments", y=0.995)
fig.tight_layout(rect=(0, 0, 0.80, 0.98))
out_path = str(snakemake.output.plot)
os.makedirs(os.path.dirname(out_path), exist_ok=True)
fig.savefig(out_path)
plt.close(fig)
