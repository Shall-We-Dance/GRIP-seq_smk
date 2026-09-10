"""Call strand-aware GRIP-seq RT-termination sites and extract sequence.

Read 2 is cDNA, so its RNA strand is opposite its alignment strand. The
cross-link nucleotide is one base upstream of the Read-2 start in RNA space.
"""

import os
import sys
from collections import defaultdict

import pysam


def revcomp(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)
log_handle = open(log_path, "w")
sys.stdout = log_handle
sys.stderr = log_handle

sample = str(snakemake.wildcards.sample)
min_mapq = int(snakemake.params.min_mapq)
min_reads = int(snakemake.params.min_reads)
neighbor_fold = float(snakemake.params.neighbor_fold)
flank = int(snakemake.params.flank)

bam = pysam.AlignmentFile(str(snakemake.input.bam), "rb")
fasta = pysam.FastaFile(str(snakemake.input.fasta))
chrom_lengths = dict(zip(fasta.references, fasta.lengths))

# Alignment strand -> counts at the first sequenced base of Read 2.
counts = defaultdict(lambda: {"+": defaultdict(int), "-": defaultdict(int)})
n_read2 = 0
for read in bam.fetch(until_eof=True):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    if not read.is_read2 or read.mapping_quality < min_mapq:
        continue
    chrom = bam.get_reference_name(read.reference_id)
    if chrom not in chrom_lengths:
        continue
    read_strand = "-" if read.is_reverse else "+"
    pos0 = read.reference_end - 1 if read.is_reverse else read.reference_start
    if 0 <= pos0 < chrom_lengths[chrom]:
        counts[chrom][read_strand][pos0] += 1
        n_read2 += 1

bam.close()
os.makedirs(os.path.dirname(str(snakemake.output.bed)), exist_ok=True)
n_sites = 0
n_sequences = 0
with open(str(snakemake.output.bed), "w") as bed, open(
    str(snakemake.output.fasta), "w"
) as seq_out:
    for chrom in fasta.references:
        for read_strand in ("+", "-"):
            strand_counts = counts.get(chrom, {}).get(read_strand, {})
            for pos0 in sorted(strand_counts):
                count = strand_counts[pos0]
                left = strand_counts.get(pos0 - 1, 0)
                right = strand_counts.get(pos0 + 1, 0)
                # The published GRIP-seq criterion is strictly >40 reads.
                if count <= min_reads or count < neighbor_fold * max(left, right):
                    continue

                # Read '+' derives from RNA '-'; Read '-' derives from RNA '+'.
                rna_strand = "-" if read_strand == "+" else "+"
                crosslink0 = pos0 + 1 if read_strand == "+" else pos0 - 1
                if not 0 <= crosslink0 < chrom_lengths[chrom]:
                    continue
                name = f"{sample}_RT_{n_sites + 1}"
                bed.write(
                    f"{chrom}\t{crosslink0}\t{crosslink0 + 1}\t{name}\t{count}"
                    f"\t{rna_strand}\t{pos0}\t{read_strand}\t{left}\t{right}\n"
                )
                n_sites += 1

                start = pos0 - flank
                end = pos0 + flank + 1
                if start < 0 or end > chrom_lengths[chrom]:
                    continue
                sequence = fasta.fetch(chrom, start, end).upper()
                if rna_strand == "-":
                    sequence = revcomp(sequence)
                if len(sequence) != 2 * flank + 1 or set(sequence) - set("ACGT"):
                    continue
                seq_out.write(
                    f">{name}|{chrom}:{pos0 + 1}|rna_strand={rna_strand}|count={count}\n"
                    f"{sequence}\n"
                )
                n_sequences += 1

fasta.close()
with open(str(snakemake.output.stats), "w") as stats:
    stats.write("sample\tusable_read2\trt_sites\tsequences\tmin_reads\tneighbor_fold\n")
    stats.write(
        f"{sample}\t{n_read2}\t{n_sites}\t{n_sequences}\t{min_reads}\t{neighbor_fold}\n"
    )

print(f"usable Read2: {n_read2}")
print(f"RT-termination sites: {n_sites}")
print(f"valid strand-oriented sequences: {n_sequences}")
log_handle.close()
