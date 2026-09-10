"""Convert GTF/GFF gene annotations to a strand-aware BED6 metagene set."""

import os
import re


def attribute(attributes, key):
    match = re.search(rf'(?:^|;\s*){re.escape(key)}[= ]+"?([^";]+)', attributes)
    return match.group(1) if match else None


gene_records = {}
transcript_records = {}
with open(str(snakemake.input.gtf)) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2].lower() not in ("gene", "transcript", "mrna"):
            continue
        chrom, _source, _feature, start, end, _score, strand, _frame, attrs = fields
        if strand not in ("+", "-"):
            continue
        feature = fields[2].lower()
        if feature == "gene":
            gene_id = attribute(attrs, "gene_id") or attribute(attrs, "ID")
        else:
            # GTF uses gene_id; GFF3 associates transcripts through Parent.
            gene_id = attribute(attrs, "gene_id") or attribute(attrs, "Parent")
            if gene_id and "," in gene_id:
                gene_id = gene_id.split(",", 1)[0]
        gene_name = attribute(attrs, "gene_name") or attribute(attrs, "Name") or gene_id
        if not gene_id:
            gene_id = f"{chrom}:{start}-{end}:{strand}"
            gene_name = gene_id
        record = (chrom, int(start) - 1, int(end), gene_name, strand)
        if feature == "gene":
            gene_records[gene_id] = record
        else:
            previous = transcript_records.get(gene_id)
            if previous is None or (record[2] - record[1]) > (previous[2] - previous[1]):
                transcript_records[gene_id] = record

# UCSC table-browser GTFs often omit explicit gene rows. In that case use the
# longest annotated transcript per gene, which avoids overweighting genes with
# many isoforms while retaining annotated TSS/TES and strand.
records = gene_records or transcript_records
if not records:
    raise ValueError("No gene/transcript/mRNA records were found in the supplied GTF/GFF.")

out_path = str(snakemake.output.bed)
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w") as out:
    for chrom, start, end, name, strand in sorted(
        records.values(), key=lambda row: (row[0], row[1], row[2], row[3])
    ):
        out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
