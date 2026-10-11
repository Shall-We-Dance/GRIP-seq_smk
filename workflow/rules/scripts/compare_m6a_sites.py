"""Compare inferred m6A sets, optional native-assembly Atlas BED6, and metagenes."""
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

if "snakemake" in globals():
    sys.path.insert(0, str(snakemake.params.script_dir))
from m6a_sites import (SiteAnnotation, bucket, consensus, load_atlas, match_count,
                      nearest_distance, read_tsv, site_key, write_tsv)
from metagene import prepare_annotation


def save_figure(fig, stem):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(str(stem)+"."+suffix, bbox_inches="tight", dpi=180)
    plt.close(fig)


def memberships(sets):
    labels = list(sets)
    union = set().union(*sets.values()) if sets else set()
    return labels, Counter(tuple(int(key in sets[label]) for label in labels) for key in union)


def upset(sets, stem, title, max_patterns=32):
    labels, patterns = memberships(sets)
    ordered = sorted(patterns.items(), key=lambda item: (-item[1], item[0]))
    write_tsv(str(stem)+".tsv", [dict(zip(labels, pattern), count=count) for pattern, count in ordered], labels+["count"])
    display = ordered[:max_patterns]
    omitted = sum(count for _, count in ordered[max_patterns:])
    fig = plt.figure(figsize=(max(7, len(display)*0.34), 4.2+len(labels)*0.23))
    grid = fig.add_gridspec(2, 1, height_ratios=[2, max(1, len(labels)*0.30)], hspace=0.06)
    top, bottom = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    if display:
        x = np.arange(len(display))
        bars = top.bar(x, [count for _, count in display], color="#276878")
        top.bar_label(bars, fontsize=7, rotation=90, padding=3)
        top.margins(y=0.3)
        for column, (pattern, _) in enumerate(display):
            bottom.scatter([column]*len(labels), range(len(labels)), color="#dddddd", s=20)
            active = np.flatnonzero(pattern)
            bottom.scatter([column]*len(active), active, color="#263d49", s=25)
            if len(active) > 1:
                bottom.plot([column, column], [active.min(), active.max()], color="#263d49")
        top.set_xlim(-0.7, len(display)-0.3)
        bottom.set_xlim(top.get_xlim())
    else:
        top.text(0.5, 0.5, "No sites / no estimable sets", ha="center", transform=top.transAxes)
    top.set_ylabel("Sites in exclusive intersection")
    top.set_xticks([])
    bottom.set_yticks(range(len(labels)), [f"{label} (n={len(sets[label]):,})" for label in labels])
    bottom.invert_yaxis()
    bottom.set_xticks([])
    bottom.set_xlabel(f"Each column is an exact membership pattern; {omitted:,} sites in omitted patterns")
    fig.suptitle(title)
    save_figure(fig, stem)


def run(smk):
    settings = dict(smk.params.settings)
    root = Path(smk.output.report)
    root.mkdir(parents=True, exist_ok=True)
    sample_names = list(smk.params.samples)
    libraries = {sample: {site_key(row) for row in read_tsv(path)} for sample, path in zip(sample_names, smk.input.calls)}
    groups = settings["groups"]
    exclude = set(settings.get("exclude_samples", []))
    minimum = int(settings.get("min_library_support", 2))
    definitions, schemes = [], {}
    for scheme, support in (("union", 1), ("replicated", minimum)):
        sets = {}
        for group, members in groups.items():
            included = [sample for sample in members if sample not in exclude]
            estimable = len(included) >= support
            sites = consensus([libraries[sample] for sample in included], support) if estimable else set()
            definitions.append(dict(scheme=scheme, group=group, minimum_support=support,
                                    included_libraries=";".join(included), excluded_libraries=";".join(s for s in members if s in exclude),
                                    status="estimable" if estimable else "not_estimable_insufficient_libraries",
                                    n_sites=len(sites) if estimable else None))
            if estimable:
                sets[group] = sites
            rows = []
            for chrom, pos, strand in sorted(sites):
                supporting = [s for s in included if (chrom, pos, strand) in libraries[s]]
                rows.append(dict(chrom=chrom, pos0=pos, pos1=pos+1, rna_strand=strand,
                                 support_libraries=";".join(supporting), n_support_libraries=len(supporting)))
            write_tsv(root/f"{scheme}.{group}.sites.tsv", rows,
                      ["chrom", "pos0", "pos1", "rna_strand", "support_libraries", "n_support_libraries"])
        schemes[scheme] = sets
        unavailable = [row["group"] for row in definitions if row["scheme"] == scheme and row["status"] != "estimable"]
        note = f"; not estimable: {', '.join(unavailable)}" if unavailable else ""
        upset(sets, root/f"UpSet.{scheme}", f"Candidate m6A: {scheme}{note}")
    write_tsv(root/"group_definitions.tsv", definitions, ["scheme", "group", "minimum_support", "included_libraries",
                                                        "excluded_libraries", "status", "n_sites"])
    upset(libraries, root/"UpSet.libraries", "Candidate m6A: individual libraries")
    all_schemes = {"libraries": libraries, **schemes}
    overlaps = []
    for scheme, sets in all_schemes.items():
        for a, b in combinations(sets, 2):
            for tolerance in (0, 1, 5):
                shared = match_count(sets[a], sets[b], tolerance)
                union = len(sets[a])+len(sets[b])-shared
                overlaps.append(dict(scheme=scheme, set_a=a, set_b=b, tolerance_nt=tolerance,
                                     n_a=len(sets[a]), n_b=len(sets[b]), shared_one_to_one=shared,
                                     a_only=len(sets[a])-shared, b_only=len(sets[b])-shared,
                                     fraction_a_in_b=shared/len(sets[a]) if sets[a] else None,
                                     fraction_b_in_a=shared/len(sets[b]) if sets[b] else None,
                                     jaccard=shared/union if union else None))
    write_tsv(root/"pairwise_overlap.tsv", overlaps, ["scheme", "set_a", "set_b", "tolerance_nt", "n_a", "n_b",
              "shared_one_to_one", "a_only", "b_only", "fraction_a_in_b", "fraction_b_in_a", "jaccard"])

    atlas, atlas_meta, annotation = None, None, None
    with pysam.FastaFile(str(smk.input.fasta)) as fasta:
        if smk.input.atlas:
            path = Path(smk.input.atlas[0])
            atlas = load_atlas(path, fasta)
            atlas_meta = dict(label=settings.get("atlas_label", "m6A-Atlas"), sites=len(atlas),
                              source=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                              coordinates="native-reference BED6, 0-based single A; RNA strand; all bases validated")
        if smk.input.gtf:
            models, selection, summary = prepare_annotation(str(smk.input.gtf[0]), available_chrom_sizes=dict(zip(fasta.references, fasta.lengths)))
            annotation = SiteAnnotation(models)
            (root/"transcript_selection.json").write_text(json.dumps(dict(summary=summary, selection=selection), indent=2)+"\n")
    coverage, annotations, profiles = [], [], []
    atlas_index = bucket(atlas) if atlas is not None else None
    for scheme, sets in all_schemes.items():
        for label, sites in sets.items():
            distances = {key: nearest_distance(key, atlas_index) for key in sites} if atlas is not None else {}
            if atlas is not None:
                bins = Counter("exact" if d == 0 else "distance1_only" if d == 1 else "distance2to5_only"
                               if d is not None and d <= 5 else "unmatched_within5" for d in distances.values())
                coverage.append(dict(scheme=scheme, label=label, n_sites=len(sites), **{name: bins[name] for name in
                                     ("exact", "distance1_only", "distance2to5_only", "unmatched_within5")}))
            for chrom, position, strand in sorted(sites):
                key = chrom, position, strand
                d = distances.get(key)
                assignments = annotation.annotate(key) if annotation else []
                annotations.append(dict(scheme=scheme, label=label, chrom=chrom, pos0=position, pos1=position+1,
                                        rna_strand=strand, genes=";".join(a[0] for a in assignments),
                                        transcripts=";".join(a[1] for a in assignments),
                                        Atlas_exact=int(d == 0) if atlas is not None else None,
                                        Atlas_within5=int(d is not None and d <= 5) if atlas is not None else None,
                                        nearest_Atlas_distance=d))
            if annotation:
                subsets = {"all": sites}
                if atlas is not None:
                    subsets.update(Atlas_exact=sites & atlas,
                                   Atlas_unmatched_within5={key for key, d in distances.items() if d is None or d > 5})
                nbins = int(settings.get("metagene_bins", 30))
                for subset, members in subsets.items():
                    counts, density, assigned = annotation.profile(members, nbins)
                    profiles += [dict(scheme=scheme, label=label, subset=subset, x=(i+0.5)/nbins,
                                      weighted_sites=count, total_sites=len(members), assigned_sites=assigned, density=density[i])
                                 for i, count in enumerate(counts)]
    annotation_fields = ["scheme", "label", "chrom", "pos0", "pos1", "rna_strand", "genes", "transcripts",
                         "Atlas_exact", "Atlas_within5", "nearest_Atlas_distance"]
    write_tsv(root/"site_annotations.tsv", annotations, annotation_fields)
    if atlas is not None:
        write_tsv(root/"Atlas_unmatched_exact.tsv", [row for row in annotations if not row["Atlas_exact"]], annotation_fields)
        write_tsv(root/"Atlas_unmatched_within5.tsv", [row for row in annotations if not row["Atlas_within5"]], annotation_fields)
    write_tsv(root/"Atlas_coverage.tsv", coverage, ["scheme", "label", "n_sites", "exact", "distance1_only", "distance2to5_only", "unmatched_within5"])
    write_tsv(root/"site_metagene.tsv", profiles, ["scheme", "label", "subset", "x", "weighted_sites", "total_sites", "assigned_sites", "density"])

    containment = []
    for scheme, groupsets in schemes.items():
        sets = dict(groupsets)
        if atlas is not None:
            union = set().union(*sets.values()) if sets else set()
            excluded = len(atlas-union)
            upset({**sets, "Atlas": atlas & union}, root/f"UpSet.Atlas.{scheme}",
                  f"{scheme}: experimental universe; {excluded:,} Atlas-only sites excluded")
            # A separate full-universe table retains database-only membership.
            labels, patterns = memberships({**sets, "Atlas": atlas})
            write_tsv(root/f"UpSet.Atlas_full.{scheme}.tsv",
                      [dict(zip(labels, pattern), count=count) for pattern, count in sorted(patterns.items())], labels+["count"])
            sets["Atlas"] = atlas
        labels = list(sets)
        matrix = np.full((len(labels), len(labels)), np.nan)
        for i, a in enumerate(labels):
            for j, b in enumerate(labels):
                shared = len(sets[a] & sets[b])
                value = 100*shared/len(sets[a]) if sets[a] else None
                matrix[i, j] = value if value is not None else np.nan
                containment.append(dict(scheme=scheme, row=a, column=b, row_total=len(sets[a]), column_total=len(sets[b]),
                                        shared=shared, percent_of_row_in_column=value))
        if labels:
            fig, ax = plt.subplots(figsize=(max(5, len(labels)*0.9), max(4, len(labels)*0.75)))
            plot = ax.imshow(matrix, vmin=0, vmax=100, cmap="Blues")
            for i in range(len(labels)):
                for j in range(len(labels)):
                    ax.text(j, i, f"{matrix[i,j]:.1f}%" if np.isfinite(matrix[i,j]) else "NA", ha="center", va="center",
                            color="white" if matrix[i,j] > 60 else "black", fontsize=9)
            ax.set_xticks(range(len(labels)), labels, rotation=35, ha="right")
            ax.set_yticks(range(len(labels)), [f"{label} (n={len(sets[label]):,})" for label in labels])
            ax.set_title(f"{scheme}: exact containment; denominator = row")
            fig.colorbar(plot, ax=ax, label="Percent of row sites")
            save_figure(fig, root/f"containment.{scheme}")
    write_tsv(root/"directional_containment.tsv", containment,
              ["scheme", "row", "column", "row_total", "column_total", "shared", "percent_of_row_in_column"])
    if profiles:
        subsets = ["all", "Atlas_exact", "Atlas_unmatched_within5"] if atlas is not None else ["all"]
        for scheme in schemes:
            fig, axes = plt.subplots(1, len(subsets), figsize=(5*len(subsets), 3.8), squeeze=False)
            for ax, subset in zip(axes[0], subsets):
                for label in schemes[scheme]:
                    rows = [r for r in profiles if r["scheme"] == scheme and r["label"] == label and r["subset"] == subset]
                    if rows and rows[0]["assigned_sites"]:
                        ax.plot([r["x"] for r in rows], [r["density"] for r in rows], label=f"{label} (n={rows[0]['assigned_sites']:,})")
                for x in (1, 2):
                    ax.axvline(x, color="grey", linestyle="--", linewidth=0.7)
                ax.set(xlim=(0, 3), ylim=(0, None), xticks=[0.5, 1.5, 2.5], xticklabels=["5′ UTR", "CDS", "3′ UTR"],
                       ylabel="Site probability density", title=f"{scheme}: {subset}")
                if ax.get_legend_handles_labels()[0]:
                    ax.legend(fontsize=8)
            save_figure(fig, root/f"site_metagene.{scheme}")
    fig, ax = plt.subplots(figsize=(7, 4))
    for sample, path in zip(sample_names, smk.input.offsets):
        rows = read_tsv(path)
        if rows and int(rows[0]["eligible_RT_sites"]):
            ax.plot([int(r["crosslink_minus_A_RNA_nt"]) for r in rows], [float(r["fraction"]) for r in rows], label=sample)
    ax.set(xlabel="Crosslink position relative to DRACH A (nt; RNA orientation)", ylabel="Fraction of unselected RT sites", ylim=(0, None))
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)
    save_figure(fig, root/"DRACH_unselected_RT_positions")
    report = dict(parameters=settings, library_counts={s: len(v) for s, v in libraries.items()}, group_definitions=definitions,
                  atlas=atlas_meta, overlap="exact same-strand adenosines; +/-1 and +/-5 use one-to-one maximum-cardinality matching",
                  atlas_coverage="same-strand any-hit proximity; unmatched is not proof of novelty",
                  metagene="unique sites with fractional multiple-model assignment; area=1; not read weighted",
                  library_stats=[json.loads(Path(p).read_text()) for p in smk.input.stats])
    (root/"report.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    run(snakemake)
