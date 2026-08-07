# workflow/scripts/average_bigwigs.py
import argparse
import os
import sys
from types import SimpleNamespace
import numpy as np
import pyBigWig

BLOCK = 1_000_000


def resolve_args():
    if "snakemake" not in globals():
        ap = argparse.ArgumentParser()
        ap.add_argument("--bws", nargs="+", required=True)
        ap.add_argument("--out-bw", required=True)
        return ap.parse_args(), None
    return (
        SimpleNamespace(bws=list(snakemake.input), out_bw=snakemake.output.bw),
        snakemake.log[0] if snakemake.log else None,
    )


def main():
    args, log_path = resolve_args()
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        log_fh = open(log_path, "w")
        sys.stdout = log_fh
        sys.stderr = log_fh

    readers = [pyBigWig.open(p) for p in args.bws]
    chroms = list(readers[0].chroms().items())
    n = len(readers)
    if n == 0:
        raise RuntimeError("No input BigWig files provided.")

    os.makedirs(os.path.dirname(args.out_bw), exist_ok=True)
    out = pyBigWig.open(args.out_bw, "w")
    out.addHeader(chroms)

    for chrom, length in chroms:
        for start in range(0, length, BLOCK):
            end = min(start + BLOCK, length)
            size = end - start
            mean = np.zeros(size, dtype=np.float64)
            for r in readers:
                vals = r.values(chrom, start, end, numpy=True)
                if vals is None:
                    continue
                mean += np.nan_to_num(vals, nan=0.0)
            mean /= n

            # write only nonzero bins; bigwigCompare treats uncovered bins as 0
            nonzero = np.flatnonzero(mean)
            if nonzero.size == 0:
                continue
            pos = nonzero + start
            out.addEntries(
                [chrom] * nonzero.size,
                pos.tolist(),
                ends=(pos + 1).tolist(),
                values=mean[nonzero].tolist(),
            )

    out.close()
    for r in readers:
        r.close()
    if log_path:
        log_fh.close()


if __name__ == "__main__":
    main()
