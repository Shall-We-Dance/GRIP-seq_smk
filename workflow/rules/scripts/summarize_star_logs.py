"""Summarize the STAR metrics most relevant to short-read retention."""

import os


def parse_star_log(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = (part.strip() for part in line.split("|", 1))
            values[key] = value.rstrip("%")
    return values


columns = [
    ("input_reads", "Number of input reads"),
    ("unique_reads", "Uniquely mapped reads number"),
    ("unique_pct", "Uniquely mapped reads %"),
    ("too_many_loci_pct", "% of reads mapped to too many loci"),
    ("too_many_mismatches_pct", "% of reads unmapped: too many mismatches"),
    ("too_short_pct", "% of reads unmapped: too short"),
]

out_path = str(snakemake.output[0])
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w") as out:
    out.write("sample\t" + "\t".join(name for name, _key in columns) + "\n")
    for path in snakemake.input:
        sample = os.path.basename(os.path.dirname(str(path)))
        metrics = parse_star_log(str(path))
        out.write(
            sample
            + "\t"
            + "\t".join(metrics.get(key, "NA") for _name, key in columns)
            + "\n"
        )
