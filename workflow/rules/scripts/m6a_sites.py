"""Motif-supported candidate m6A inference; no database-based site selection."""
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import sys

if "snakemake" in globals():
    sys.path.insert(0, str(snakemake.params.script_dir))
from call_rt_termination_sites import mapped_crosslinks, stream_selected_sites
from motif_utils import revcomp, text_open

DRACH = re.compile("[AGT][AG]AC[ACT]")
SITE_FIELDS = ["chrom", "pos0", "pos1", "rna_strand", "sample", "crosslink_pos0",
               "r2_positions0", "RT_peak_support_reads", "left_neighbor_reads", "right_neighbor_reads",
               "crosslink_reads", "Input_crosslink_reads", "IP_Input_fold", "RNA_5mer",
               "A_at_offset", "DRACH_at_offset", "passes_RT_Input", "passes_m6A"]


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def site_key(row):
    return row["chrom"], int(row["pos0"]), row["rna_strand"]


def sequence_at(fasta, key, flank):
    chrom, position, strand = key
    if chrom not in fasta.references or position < flank or position + flank >= fasta.get_reference_length(chrom):
        return None
    sequence = fasta.fetch(chrom, position-flank, position+flank+1).upper()
    return revcomp(sequence) if strand == "-" else sequence


def counts_at(path, keys):
    result = {}
    with text_open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, _, count, strand = line.rstrip().split("\t")[:6]
            key = chrom, int(start), strand
            if key in keys:
                if key in result or int(end) != key[1]+1 or int(count) < 1:
                    raise ValueError("Invalid or duplicate counted crosslink record")
                result[key] = int(count)
    return result


def enrichment(ip, control, ip_depth, control_depth, pseudocount=0.5):
    if min(ip_depth, control_depth) <= 0:
        return None
    return ((ip+pseudocount)/ip_depth) / ((control+pseudocount)/control_depth)


def infer_sites(sites, ip_counts, input_counts, fasta, sample, ip_depth, input_depth, settings):
    rows, offsets, eligible = [], Counter(), 0
    offset = int(settings.get("adenosine_offset", 3))
    window = int(settings.get("offset_window", 25))
    for key, support in sorted(sites.items()):
        chrom, crosslink, strand = key
        direction = 1 if strand == "+" else -1
        adenosine = chrom, crosslink+direction*offset, strand
        five = sequence_at(fasta, adenosine, 2) or ""
        is_a = bool(five and five[2] == "A")
        drach = bool(DRACH.fullmatch(five))
        ip, control = ip_counts.get(key, 0), input_counts.get(key, 0)
        if support["count"] > ip:
            raise ValueError(f"RT support exceeds crosslink count at {key}")
        fold = enrichment(ip, control, ip_depth, input_depth, float(settings.get("pseudocount", 0.5)))
        passes = (support["count"] > int(settings.get("min_reads", 10)) and
                  fold is not None and fold >= float(settings.get("min_input_fold", 2)))
        rows.append(dict(chrom=chrom, pos0=adenosine[1], pos1=adenosine[1]+1, rna_strand=strand,
                         sample=sample, crosslink_pos0=crosslink,
                         r2_positions0=",".join(map(str, sorted(support["r2_positions"]))),
                         RT_peak_support_reads=support["count"], left_neighbor_reads=support["left"],
                         right_neighbor_reads=support["right"], crosslink_reads=ip, Input_crosslink_reads=control,
                         IP_Input_fold=fold, RNA_5mer=five, A_at_offset=int(is_a), DRACH_at_offset=int(drach),
                         passes_RT_Input=int(passes),
                         passes_m6A=int(passes and (drach if settings.get("require_DRACH", True) else is_a))))
        # All RT loci enter this diagnostic before Input and motif selection.
        sequence = sequence_at(fasta, key, window+2)
        if sequence and not set(sequence)-set("ACGT"):
            eligible += 1
            for i in range(len(sequence)-4):
                if DRACH.fullmatch(sequence[i:i+5]):
                    offsets[window-i] += 1  # crosslink minus central A along RNA
    return rows, offsets, eligible


def bucket(sites):
    result = defaultdict(list)
    for chrom, position, strand in sites:
        result[chrom, strand].append(position)
    return {key: sorted(set(values)) for key, values in result.items()}


def nearest_distance(key, index):
    chrom, position, strand = key
    values = index.get((chrom, strand), [])
    i = bisect_left(values, position)
    near = values[max(0, i-1):i+1]
    return min((abs(x-position) for x in near), default=None)


def match_count(first, second, tolerance=0):
    if not tolerance:
        return len(first & second)
    a, b, matched = bucket(first), bucket(second), 0
    for key, x in a.items():
        y, i, j = b.get(key, []), 0, 0
        while i < len(x) and j < len(y):
            if x[i] < y[j]-tolerance:
                i += 1
            elif y[j] < x[i]-tolerance:
                j += 1
            else:
                matched += 1
                i += 1
                j += 1
    return matched


def consensus(sets, minimum):
    return {key for key, count in Counter(key for sites in sets for key in sites).items() if count >= minimum}


def load_atlas(path, fasta):
    """Read native-assembly BED6 single-A coordinates, never guess a coordinate system."""
    sites = set()
    with text_open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 6:
                raise ValueError("Atlas requires BED6: chromosome, start0, end, name, score, RNA strand")
            chrom, start, end, _, _, strand = fields[:6]
            key = chrom, int(start), strand
            if strand not in {"+", "-"} or key[1] < 0 or int(end) != key[1]+1 or sequence_at(fasta, key, 0) != "A":
                raise ValueError(f"Atlas record is not a 1-nt RNA-strand adenosine in this reference: {chrom}:{start}")
            sites.add(key)
    return sites


class SiteAnnotation:
    def __init__(self, models):
        self.models, self.index = models, defaultdict(set)
        for i, model in enumerate(models):
            for start, end in model["exons"]:
                for block in range(start//10000, (end-1)//10000+1):
                    self.index[model["chrom"], model["strand"], block].add(i)

    def annotate(self, key):
        chrom, position, strand = key
        found = []
        for i in sorted(self.index.get((chrom, strand, position//10000), [])):
            model = self.models[i]
            exons = model["exons"] if strand == "+" else model["exons"][::-1]
            offset = 0
            for start, end in exons:
                if start <= position < end:
                    offset += position-start if strand == "+" else end-1-position
                    u, c, v = (model["lengths"][name] for name in ("utr5", "cds", "utr3"))
                    scaled = offset/u if offset < u else 1+(offset-u)/c if offset < u+c else 2+(offset-u-c)/v
                    found.append((model["gene_id"], model["transcript_id"], scaled))
                    break
                offset += end-start
        return [(gene, transcript, scaled, 1/len(found)) for gene, transcript, scaled in found]

    def profile(self, sites, bins=30):
        counts, assigned = [0.0]*(3*bins), 0
        for key in sites:
            annotations = self.annotate(key)
            assigned += bool(annotations)
            for _, _, position, weight in annotations:
                counts[min(3*bins-1, int(position*bins))] += weight
        density = [count/assigned*bins if assigned else None for count in counts]
        return counts, density, assigned


def run(smk):
    import pysam
    settings = dict(smk.params.settings)
    selected, counts_stats = stream_selected_sites(smk.input.counts, int(settings.get("min_reads", 10)),
                                                  float(settings.get("neighbor_fold", 1.5)))
    sites, mapping_stats = mapped_crosslinks(smk.input.mapping, selected, int(smk.params.crosslink_shift))
    ip_counts = counts_at(smk.input.crosslinks, sites)
    input_counts = counts_at(smk.input.control, sites)
    depths = [json.loads(Path(path).read_text())["accepted_read2"] for path in (smk.input.ip_stats, smk.input.control_stats)]
    with pysam.FastaFile(str(smk.input.fasta)) as fasta:
        rows, offsets, eligible = infer_sites(sites, ip_counts, input_counts, fasta, str(smk.wildcards.sample), *depths, settings)
    candidates = [row for row in rows if row["passes_m6A"]]
    write_tsv(smk.output.all_rt, rows, SITE_FIELDS)
    write_tsv(smk.output.candidates, candidates, SITE_FIELDS)
    with open(smk.output.bed, "w") as handle:
        for i, row in enumerate(candidates, 1):
            handle.write(f"{row['chrom']}\t{row['pos0']}\t{row['pos0']+1}\t{row['sample']}_{i}\t0\t{row['rna_strand']}\n")
    window = int(settings.get("offset_window", 25))
    write_tsv(smk.output.offsets, [dict(crosslink_minus_A_RNA_nt=x, site_anchor_pairs=offsets[x],
              eligible_RT_sites=eligible, fraction=offsets[x]/eligible if eligible else None) for x in range(-window, window+1)],
              ["crosslink_minus_A_RNA_nt", "site_anchor_pairs", "eligible_RT_sites", "fraction"])
    report = dict(sample=str(smk.wildcards.sample), control=str(smk.params.control), parameters=settings,
                  variant=str(smk.params.variant), RT_sites=len(sites), candidates=len(candidates),
                  IP_accepted_R2=depths[0], Input_accepted_R2=depths[1], **dict(counts_stats), **dict(mapping_stats),
                  status="ok" if min(depths) > 0 else "not_estimable_zero_IP_or_Input_depth",
                  normalization="counts from selected variant; denominators are accepted R2 before blacklist",
                  site_definition=("DRACH-supported" if settings.get("require_DRACH", True) else "A-only") + " inferred m6A; not chemical confirmation",
                  sequence_context="strand-oriented contiguous genomic coordinates; not splice-junction reconstruction")
    Path(smk.output.stats).write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    run(snakemake)
