# workflow/scripts/motif_anchor.py
# Anchor R2-first-base 1-nt sites to an RNA-motif (default DRACH, m6A = 3rd base A).
# Per site:
#   1) shift to the cross-link/m6A site (RNA 5'-direction; fwd: -1, rev: +1)
#   2) optional motif anchoring: exact match, or nearest match within a window
#   3) output a corrected 1-nt BED (BED6 + motif sequence column)
import argparse
import os
import re
import sys
from types import SimpleNamespace
import pysam

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def revcomp(seq):
    return "".join(COMP.get(b, "N") for b in reversed(seq.upper()))


def load_fai(fai_path):
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chroms.append((fields[0], int(fields[1])))
    return chroms


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed-in", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--motif", default="DRACH")
    ap.add_argument("--motif-pos", type=int, default=3)
    ap.add_argument("--site-shift", type=int, default=1)
    ap.add_argument("--search-window", type=int, default=0)
    ap.add_argument("--keep-unmatched", action="store_true")
    return ap.parse_args()


def resolve_args():
    if "snakemake" not in globals():
        return parse_args(), None
    log_path = snakemake.log[0] if snakemake.log else None
    args = SimpleNamespace(
        bed_in=snakemake.input.bed,
        fasta=snakemake.input.fasta,
        fai=snakemake.input.fai,
        out_bed=snakemake.output.bed,
        motif=str(snakemake.params.motif),
        motif_pos=int(snakemake.params.motif_pos),
        site_shift=int(snakemake.params.site_shift),
        search_window=int(snakemake.params.search_window),
        keep_unmatched=bool(snakemake.params.keep_unmatched),
    )
    return args, log_path


def main():
    args, log_path = resolve_args()
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        log_fh = open(log_path, "w")
        sys.stdout = log_fh
        sys.stderr = log_fh

    motif = args.motif.upper()
    motif_len = len(motif)
    if not 1 <= args.motif_pos <= motif_len:
        raise ValueError(f"motif_pos ({args.motif_pos}) must be within 1..{motif_len} of motif '{motif}'")
    motif_re = re.compile("".join(IUPAC.get(b, "N") for b in motif))

    half5 = args.motif_pos - 1          # nt upstream (RNA 5' side) of the site within the motif
    half3 = motif_len - args.motif_pos  # nt downstream (RNA 3' side)

    chrom_len = dict(load_fai(args.fai))
    fasta = pysam.FastaFile(args.fasta)

    def fetch_around(chrom, site, strand):
        """
        Extract motif-length RNA-orientation sequence centered on `site`.
        Read2 is cDNA (reverse complement of the RNA); RNA strand is the
        opposite of the read strand: read '+' => RNA '-', read '-' => RNA '+'.
        """
        L = chrom_len.get(chrom)
        if L is None:
            return None
        if strand == "-":
            begin = site - half5
            end = site + half3 + 1
        else:
            begin = site - half3
            end = site + half5 + 1
        if begin < 0 or end > L:
            return None
        seq = fasta.fetch(chrom, begin, end).upper()
        if strand == "+":
            seq = revcomp(seq)
        return seq if len(seq) == motif_len else None

    def genomic_offset(strand, offset_rna):
        return offset_rna if strand == "+" else -offset_rna

    n_in = 0
    n_anchored = 0
    n_unmatched_kept = 0
    n_dropped = 0

    os.makedirs(os.path.dirname(args.out_bed), exist_ok=True)
    with open(args.bed_in) as fin, open(args.out_bed, "w") as fout:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, start0, end, name, score, strand = fields[:6]
            if chrom not in chrom_len:
                continue
            n_in += 1
            pos = int(start0)
            # Read2 is cDNA (RNA is the opposite strand): the m6A/cross-link site is
            # +site_shift for '+' reads (RNA '-'), -site_shift for '-' reads (RNA '+').
            site = pos + (args.site_shift if strand == "+" else -args.site_shift)

            seq = fetch_around(chrom, site, strand)
            found = seq is not None and motif_re.fullmatch(seq)
            anchored = site
            if not found and args.search_window > 0:
                for d in range(1, args.search_window + 1):
                    for off in (d, -d):
                        cand = site + genomic_offset(strand, off)
                        cseq = fetch_around(chrom, cand, strand)
                        if cseq is not None and motif_re.fullmatch(cseq):
                            anchored = cand
                            seq = cseq
                            found = True
                            break
                    if found:
                        break

            if not found:
                if args.keep_unmatched:
                    n_unmatched_kept += 1
                    fout.write(
                        f"{chrom}\t{site}\t{site + 1}\t{name}\t{score}\t{strand}\tNA\n"
                    )
                else:
                    n_dropped += 1
                continue

            n_anchored += 1
            fout.write(
                f"{chrom}\t{anchored}\t{anchored + 1}\t{name}\t{score}\t{strand}\t{seq}\n"
            )

    fasta.close()

    print(f"[motif_anchor] motif={motif} motif_pos={args.motif_pos} "
          f"site_shift={args.site_shift} search_window={args.search_window} "
          f"keep_unmatched={args.keep_unmatched}")
    print(f"[motif_anchor] sites in: {n_in}")
    print(f"[motif_anchor] anchored: {n_anchored}")
    print(f"[motif_anchor] unmatched kept: {n_unmatched_kept}")
    print(f"[motif_anchor] unmatched dropped: {n_dropped}")

    if log_path:
        log_fh.close()


if __name__ == "__main__":
    main()
