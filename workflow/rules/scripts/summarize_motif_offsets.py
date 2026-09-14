"""Summarize every sample motif, with group plots restricted to predefined motifs.

De novo consensuses remain sample-specific, even when their MEME names match;
the summary retains their actual IUPAC definition, anchor and validation results.

Plots show equal-sample means of per-sample site fractions. SD bands describe
variation among samples with effective sites, not confidence in a pooled-site
estimate. Zero-site samples are NA and excluded from means; valid samples with
zero motif matches contribute observed zero fractions. Existing sample-level
p/q values are copied for inspection and are never recombined across samples.
"""

import csv
import json
import math
from pathlib import Path


SUMMARY_FIELDS = [
    "sample", "comparison_roles", "motif", "iupac", "anchor_pos_1based", "source", "status",
    "effective_site_pairs", "supporting_reads", "peak_frequency_offsets_nt", "peak_site_fraction",
    "qvalues_at_frequency_peaks", "peak_excess_over_background_offsets_nt", "minimum_qvalue_bh_sample",
    "significant_enrichment_offsets_q_lt_0p05", "offset_zero_foreground_fraction",
    "offset_zero_background_fraction", "offset_zero_qvalue_bh_sample", "qvalue_scope",
]
CURVE_FIELDS = [
    "group", "role", "motif", "iupac", "anchor_pos_1based",
    "motif_anchor_offset_from_crosslink_rna", "population", "mean_site_fraction", "sd_samples",
    "n_samples_effective", "n_samples_configured",
]


def numeric(value):
    if value is None or value in ("", "NA", "nan", "NaN"):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def write_rows(path, columns, rows):
    with open(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if value is None else value for key, value in row.items()})


def read_results(sample_names, tables, statistics):
    """Validate matching definitions/coordinates before comparing any samples."""
    if len(sample_names) != len(tables) or len(sample_names) != len(statistics):
        raise ValueError("Samples, motif tables, and statistics must have matching lengths")
    profiles, reports, definitions, offsets = {}, {}, {}, {}
    for sample, table, stat_path in zip(sample_names, tables, statistics):
        with open(stat_path) as handle:
            report = json.load(handle)
        if report.get("sample") != sample:
            raise ValueError(f"Expected motif statistics for {sample}, got {report.get('sample')}")
        reports[sample] = report
        motif_stats = {motif["name"]: motif for motif in report["motifs"]}
        if len(motif_stats) != len(report["motifs"]):
            raise ValueError(f"{sample}: duplicate motif names in statistics")
        with open(table) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row["motif"]
                if name not in motif_stats:
                    raise ValueError(f"{sample}: table motif {name} absent from statistics")
                definition = (row["iupac"], int(row["anchor_pos_1based"]))
                summary = motif_stats[name]
                if row["source"] != summary["source"]:
                    raise ValueError(f"{sample}: motif {name} source differs between table/statistics")
                if row["source"] == "predefined":
                    if name in definitions and definitions[name] != definition:
                        raise ValueError(f"Predefined motif {name} differs between samples")
                    definitions[name] = definition
                if definition != (summary["sequence"], int(summary["anchor_pos"])):
                    raise ValueError(f"{sample}: motif {name} definition differs between table/statistics")
                n_sites = int(row["effective_site_pairs"])
                if n_sites != int(summary["effective_site_pairs"]):
                    raise ValueError(f"{sample}: inconsistent effective denominator for {name}")
                if row.get("coordinate_population", "all_eligible_crosslinks") != "all_eligible_crosslinks":
                    raise ValueError("Motif group summaries require crosslink-population offsets.tsv")
                offset = int(row["motif_anchor_offset_from_crosslink_rna"])
                key = (sample, name)
                if offset in profiles.setdefault(key, {}):
                    raise ValueError(f"Duplicate motif offset in {sample}: {name} {offset}")
                clean = dict(row)
                clean["offset"] = offset
                for column in ("foreground_site_fraction", "background_site_fraction", "qvalue_bh_sample"):
                    clean[column] = numeric(row[column]) if n_sites > 0 else None
                if n_sites > 0 and any(clean[c] is None for c in ("foreground_site_fraction", "background_site_fraction")):
                    raise ValueError(f"{sample}: effective motif sites require finite foreground/background fractions")
                profiles[key][offset] = clean
        for name in motif_stats:
            key = (sample, name)
            if key not in profiles:
                raise ValueError(f"{sample}: motif {name} has statistics but no offset rows")
            current = sorted(profiles[key])
            if motif_stats[name]["source"] == "predefined":
                if name in offsets and current != offsets[name]:
                    raise ValueError(f"Offset grid for {name} differs between samples")
                offsets[name] = current
    return profiles, reports, definitions, offsets


def comparison_groups(sample_names, configured):
    groups = {
        name: {role: [sample for sample in members if sample in sample_names]
               for role, members in roles.items() if role in ("IP", "Input")}
        for name, roles in configured.items()
    }
    groups = {name: roles for name, roles in groups.items() if any(roles.values())}
    assigned = {sample for roles in groups.values() for members in roles.values() for sample in members}
    ungrouped = [sample for sample in sample_names if sample not in assigned]
    if ungrouped:
        groups["Individual samples"] = {sample: [sample] for sample in ungrouped}
    return groups


def sample_summary_rows(sample_names, reports, profiles, groups):
    rows = []
    for sample in sample_names:
        membership = [f"{group}:{role}" for group, roles in groups.items()
                      for role, members in roles.items() if sample in members]
        report = reports[sample]
        for motif in report["motifs"]:
            name, n = motif["name"], int(motif["effective_site_pairs"])
            table = profiles[(sample, name)]
            peaks = motif.get("peak_frequency_offsets_nt", []) if n else []
            peak_q = {str(offset): table[offset]["qvalue_bh_sample"] for offset in peaks}
            fractions = [row["foreground_site_fraction"] for row in table.values()
                         if row["foreground_site_fraction"] is not None]
            qvalues = [row["qvalue_bh_sample"] for row in table.values() if row["qvalue_bh_sample"] is not None]
            center = table.get(0, {})
            rows.append(dict(zip(SUMMARY_FIELDS, [
                sample, "|".join(membership), name, motif["sequence"], motif["anchor_pos"], motif["source"],
                motif["status"], n, motif.get("supporting_reads", 0), json.dumps(peaks),
                max(fractions) if fractions else None,
                json.dumps(peak_q, sort_keys=True) if n else None,
                json.dumps(motif.get("peak_excess_over_background_offsets_nt", []) if n else []),
                min(qvalues) if qvalues else None,
                json.dumps(motif.get("significant_enrichment_offsets_q_lt_0p05", []) if n else []),
                center.get("foreground_site_fraction"), center.get("background_site_fraction"),
                center.get("qvalue_bh_sample"), report.get("multiple_testing", "sample-level BH; copied from input"),
            ])))
    return rows


def run(smk):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import numpy as np

    sample_names = list(smk.params.samples)
    profiles, reports, definitions, offsets = read_results(
        sample_names, list(smk.input.tables), list(smk.input.statistics)
    )
    groups = comparison_groups(sample_names, dict(smk.params.groups))
    for path in smk.output:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    summary_rows = sample_summary_rows(sample_names, reports, profiles, groups)
    write_rows(smk.output.summary, SUMMARY_FIELDS, summary_rows)
    curves = []
    fig, axes = plt.subplots(max(1, len(definitions)), max(1, len(groups)),
                             figsize=(5.6 * max(1, len(groups)), 3.9 * max(1, len(definitions))),
                             squeeze=False)
    if not definitions or not groups:
        axes[0, 0].text(0.5, 0.5, "No predefined motif sample profiles available", ha="center", va="center")
        axes[0, 0].set_axis_off()
    for motif_index, (motif, definition) in enumerate(definitions.items()):
        for group_index, (group, roles) in enumerate(groups.items()):
            ax = axes[motif_index, group_index]
            missing_samples = []
            x = np.asarray(offsets[motif])
            for role_index, (role, members) in enumerate(roles.items()):
                if not members:
                    continue
                color = {"IP": "#bb463c", "Input": "#316d9d"}.get(role, plt.cm.tab20(role_index % 20))
                for population, field in (("foreground", "foreground_site_fraction"),
                                           ("local_background", "background_site_fraction")):
                    matrix = np.asarray([
                        [profiles.get((sample, motif), {}).get(offset, {}).get(field) for offset in x]
                        for sample in members
                    ], dtype=float)
                    valid = np.isfinite(matrix)
                    n = valid.sum(axis=0)
                    mean = np.divide(np.nansum(matrix, axis=0), n,
                                     out=np.full(len(x), np.nan), where=n > 0)
                    delta = np.where(valid, matrix - mean, 0)
                    sd = np.sqrt(np.divide((delta * delta).sum(axis=0), n - 1,
                                           out=np.full(len(x), np.nan), where=n > 1))
                    n_effective = int(valid.any(axis=1).sum())
                    for i, offset in enumerate(x):
                        curves.append(dict(zip(CURVE_FIELDS, [
                            group, role, motif, definition[0], definition[1], int(offset), population,
                            float(mean[i]) if np.isfinite(mean[i]) else None,
                            float(sd[i]) if np.isfinite(sd[i]) else None, int(n[i]), len(members),
                        ])))
                    if population == "foreground":
                        for sample, values in zip(members, matrix):
                            if np.isfinite(values).any():
                                ax.plot(x, values, color=color, alpha=0.22, linewidth=0.7)
                            else:
                                missing_samples.append(sample)
                        if n_effective:
                            ax.plot(x, mean, color=color, linewidth=1.8,
                                    label=f"{role} ({n_effective}/{len(members)} samples)")
                        else:
                            ax.plot([], [], color=color, linewidth=1.8, label=f"{role}: NA (no effective sites)")
                        if np.any(n > 1):
                            ax.fill_between(x, np.maximum(0, mean - sd), np.minimum(1, mean + sd),
                                            color=color, alpha=0.16, linewidth=0)
                    elif n_effective:
                        ax.plot(x, mean, color=color, linestyle="--", linewidth=1.1, alpha=0.8,
                                label=f"{role} local background")
            ax.axvline(0, color="0.6", linestyle=":", linewidth=0.8)
            ax.set_title(f"{group}\n{motif}: {definition[0]} (anchor {definition[1]})", fontsize=10)
            ax.set_ylabel("Fraction of distinct crosslink sites")
            ax.set_xlabel("Motif anchor relative to crosslink (nt; RNA 3′ positive)", fontsize=9)
            ax.set_ylim(bottom=0)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(axis="y", alpha=0.17)
            ax.spines[["top", "right"]].set_visible(False)
            ax.legend(fontsize=7, frameon=False)
            if missing_samples:
                ax.text(0.01, 0.98, "NA: " + ", ".join(missing_samples), transform=ax.transAxes,
                        ha="left", va="top", fontsize=7, color="0.4", wrap=True)
    fig.suptitle("Predefined motifs across samples\nEqual-sample means; shading = sample SD (n ≥ 2)", fontsize=12)
    fig.text(0.5, 0.008, "Faint lines: individual samples. Dashed: mean matched background. "
             "Zero-site samples are NA. No pooled-site tests or combined q-values.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fig.savefig(smk.output.plot, bbox_inches="tight")
    plt.close(fig)
    write_rows(smk.output.curves, CURVE_FIELDS, curves)


if __name__ == "__main__":
    run(snakemake)
