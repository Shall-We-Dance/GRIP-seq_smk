"""Optional motif-conditioned diagnostic; not unbiased modification discovery.

Input BED6 is on the RNA strand, with already inferred cross-link coordinates.
The workflow applies no additional site shift. BED output adds motif sequence,
original cross-link coordinate, displacement along RNA and diagnostic status.
Ambiguous equidistant hits are never resolved by an arbitrary strand bias.
"""
import argparse
import gzip
import json
import os
import re
from collections import Counter
from contextlib import redirect_stdout, redirect_stderr
from types import SimpleNamespace

import pysam

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def revcomp(seq):
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def annotate(args):
    motif = args.motif.upper()
    if not motif or set(motif) - IUPAC.keys():
        raise ValueError("motif must contain valid IUPAC DNA/RNA letters")
    if not 1 <= args.motif_pos <= len(motif) or args.search_window < 0:
        raise ValueError("Invalid motif_pos or negative search_window")
    motif_re = re.compile("".join(IUPAC[b] for b in motif))
    half5, half3 = args.motif_pos - 1, len(motif) - args.motif_pos
    os.makedirs(os.path.dirname(os.path.abspath(args.out_bed)), exist_ok=True)
    stats = Counter()
    opener = gzip.open if str(args.bed_in).endswith(".gz") else open
    with pysam.FastaFile(args.fasta) as fasta, opener(args.bed_in, "rt") as fin, open(args.out_bed, "w") as fout:
        lengths = dict(zip(fasta.references, fasta.lengths))

        def fetch(chrom, site, strand):
            lo, hi = (site - half5, site + half3 + 1) if strand == "+" else (site - half3, site + half5 + 1)
            if lo < 0 or hi > lengths[chrom]:
                return None
            seq = fasta.fetch(chrom, lo, hi).upper()
            return seq if strand == "+" else revcomp(seq)

        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, name, score, strand = line.split()[:6]
            if strand not in ("+", "-") or int(end) != int(start) + 1:
                raise ValueError("Expected 1-nt BED6 with explicit RNA strand")
            stats["input_sites"] += 1
            stats["input_reads"] += int(score)
            if chrom not in lengths:
                stats["unknown_reference"] += 1
                continue
            site = int(start) - args.site_shift * (1 if strand == "+" else -1)
            if not 0 <= site < lengths[chrom]:
                stats["reference_boundary"] += 1
                continue
            hits = []
            ambiguous = False
            for distance in range(args.search_window + 1):
                for offset in ((0,) if distance == 0 else (-distance, distance)):
                    candidate = site + offset * (1 if strand == "+" else -1)
                    seq = fetch(chrom, candidate, strand)
                    if seq is not None and motif_re.fullmatch(seq):
                        hits.append((candidate, seq, offset))
                if hits:
                    ambiguous = len(hits) > 1
                    break
            status = "ambiguous" if ambiguous else "matched" if hits else "unmatched"
            stats[status + "_sites"] += 1
            stats[status + "_reads"] += int(score)
            if status == "matched":
                chosen, sequence, offset = hits[0]
                fout.write(f"{chrom}\t{chosen}\t{chosen + 1}\t{name}\t{score}\t{strand}\t{sequence}\t{site}\t{offset}\t{status}\n")
            elif args.keep_unmatched:
                fout.write(f"{chrom}\t{site}\t{site + 1}\t{name}\t{score}\t{strand}\tNA\t{site}\t0\t{status}\n")
    result = dict(stats)
    result.update(motif=motif, motif_pos=args.motif_pos, search_window=args.search_window, site_shift=args.site_shift, strand="RNA", interpretation="motif-conditioned diagnostic only; neither unbiased localization nor confirmed modification sites", sequence_context="genomic; windows may include intronic sequence near exon boundaries")
    if getattr(args, "out_stats", None):
        with open(args.out_stats, "w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    print(json.dumps(result, indent=2))
    return result


def resolve_args():
    if "snakemake" in globals():
        return SimpleNamespace(bed_in=str(snakemake.input.bed), fasta=str(snakemake.input.fasta), fai=str(snakemake.input.fai), out_bed=str(snakemake.output.bed), out_stats=snakemake.output.get("stats"), motif=str(snakemake.params.motif), motif_pos=int(snakemake.params.motif_pos), site_shift=int(snakemake.params.site_shift), search_window=int(snakemake.params.search_window), keep_unmatched=bool(snakemake.params.keep_unmatched)), str(snakemake.log[0]) if snakemake.log else None
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("bed-in", "fasta", "fai", "out-bed"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--out-stats")
    parser.add_argument("--motif", default="DRACH")
    parser.add_argument("--motif-pos", type=int, default=3)
    parser.add_argument("--site-shift", type=int, default=0, help="RNA-upstream shift; workflow input is already shifted")
    parser.add_argument("--search-window", type=int, default=0)
    parser.add_argument("--keep-unmatched", action="store_true")
    return parser.parse_args(), None


if __name__ == "__main__":
    args, log_path = resolve_args()
    if log_path:
        os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
        with open(log_path, "w") as handle, redirect_stdout(handle), redirect_stderr(handle):
            annotate(args)
    else:
        annotate(args)
