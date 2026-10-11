"""Audit the first retained paired FASTQ records without random-sampling claims."""
import argparse
import gzip
import json


def audit(r1, r2, sample, minimum_length=20, adapter_prefix="AGATCGGAAGAGC", limit=20000):
    n, adapters, minimum = 0, [0, 0], [None, None]
    with gzip.open(r1, "rt") as first, gzip.open(r2, "rt") as second:
        for _ in range(limit):
            pair = [[handle.readline().rstrip("\r\n") for _ in range(4)] for handle in (first, second)]
            if not pair[0][0] and not pair[1][0]:
                break
            ids = [row[0].split()[0].removesuffix(f"/{i+1}") if row[0] else "" for i, row in enumerate(pair)]
            if not ids[0] or ids[0] != ids[1]:
                raise ValueError("Paired FASTQ names or record counts do not match")
            for i, row in enumerate(pair):
                if not (row[0].startswith("@") and row[2].startswith("+") and
                        len(row[1]) == len(row[3]) and len(row[1]) >= minimum_length):
                    raise ValueError("Malformed or undersized retained FASTQ record")
                minimum[i] = len(row[1]) if minimum[i] is None else min(minimum[i], len(row[1]))
                adapters[i] += adapter_prefix in row[1].upper()
            n += 1
    return dict(sample=sample, sampled_pairs=n, selection=f"first {limit} retained pairs; not a random sample",
                adapter_prefix=adapter_prefix, adapter_reads_R1=adapters[0], adapter_reads_R2=adapters[1],
                minimum_lengths=minimum, minimum_length_required=minimum_length)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sample", "r1", "r2", "output"):
        parser.add_argument(name)
    parser.add_argument("--minimum-length", type=int, default=20)
    parser.add_argument("--adapter-prefix", default="AGATCGGAAGAGC")
    args = parser.parse_args()
    with open(args.output, "w") as handle:
        json.dump(audit(args.r1, args.r2, args.sample, args.minimum_length, args.adapter_prefix), handle, indent=2)
        handle.write("\n")
