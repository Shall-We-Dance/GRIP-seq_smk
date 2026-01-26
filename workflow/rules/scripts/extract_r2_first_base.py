# workflow/scripts/extract_r2_first_base.py
import argparse
import os
import sys
from types import SimpleNamespace
import pysam
import pyBigWig

def load_fai(fai_path):
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chroms.append((fields[0], int(fields[1])))
    return chroms

def load_blacklist(bed_path):
    """
    Load BED into per-chrom interval lists, then binary search membership.
    Assumes BED is 0-based half-open.
    """
    from bisect import bisect_right

    bl = {}
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, s, e = line.rstrip("\n").split("\t")[:3]
            s, e = int(s), int(e)
            bl.setdefault(chrom, []).append((s, e))

    # sort and also prepare starts for bisect
    bl_idx = {}
    for chrom, intervals in bl.items():
        intervals.sort()
        starts = [s for s, _ in intervals]
        bl_idx[chrom] = (starts, intervals)

    def is_blacklisted(chrom, pos0):
        if chrom not in bl_idx:
            return False
        starts, intervals = bl_idx[chrom]
        i = bisect_right(starts, pos0) - 1
        if i < 0:
            return False
        s, e = intervals[i]
        return (s <= pos0 < e)

    return is_blacklisted

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out-bw", required=True)
    ap.add_argument("--min-mapq", type=int, default=11)
    ap.add_argument("--blacklist-bed", default=None)
    ap.add_argument("--require-proper-pair", action="store_true")
    return ap.parse_args()

def resolve_args():
    if "snakemake" not in globals():
        return parse_args(), None

    log_path = snakemake.log[0] if snakemake.log else None
    args = SimpleNamespace(
        bam=snakemake.input.bam,
        fai=snakemake.input.fai,
        out_bw=snakemake.output.bw,
        min_mapq=int(snakemake.params.min_mapq),
        blacklist_bed=snakemake.input.get("bl"),
        require_proper_pair=snakemake.params.get("require_proper_pair", False),
    )
    return args, log_path

def main():
    args, log_path = resolve_args()
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        log_fh = open(log_path, "w")
        sys.stdout = log_fh
        sys.stderr = log_fh

    chroms = load_fai(args.fai)
    chrom_len = {c: L for c, L in chroms}

    is_blacklisted = None
    if args.blacklist_bed:
        is_blacklisted = load_blacklist(args.blacklist_bed)

    bam = pysam.AlignmentFile(args.bam, "rb")

    # count positions as sparse dict: (chrom -> {pos0: count})
    counts = {c: {} for c, _ in chroms}
    n_read2 = 0

    for r in bam.fetch(until_eof=True):
        if r.is_unmapped:
            continue
        if r.is_secondary or r.is_supplementary:
            continue
        if r.mapping_quality < args.min_mapq:
            continue
        if args.require_proper_pair and (not r.is_proper_pair):
            continue
        if not r.is_read2:
            continue

        chrom = bam.get_reference_name(r.reference_id)
        if chrom not in counts:
            continue

        # "first base pair" of Read2 in sequencing orientation:
        # if aligned forward, first base is reference_start
        # if aligned reverse, first base maps to reference_end - 1
        if r.is_reverse:
            pos0 = r.reference_end - 1
        else:
            pos0 = r.reference_start

        if pos0 < 0 or pos0 >= chrom_len.get(chrom, 0):
            continue

        if is_blacklisted and is_blacklisted(chrom, pos0):
            continue

        n_read2 += 1
        d = counts[chrom]
        d[pos0] = d.get(pos0, 0) + 1

    bam.close()

    if n_read2 == 0:
        raise RuntimeError("No usable Read2 records after filtering; cannot compute CPM.")

    scale = 1e6 / n_read2

    os.makedirs(os.path.dirname(args.out_bw), exist_ok=True)
    bw = pyBigWig.open(args.out_bw, "w")
    bw.addHeader(chroms)

    # write as 1bp intervals
    for chrom, _L in chroms:
        d = counts.get(chrom, {})
        if not d:
            continue
        positions = sorted(d.keys())
        starts = positions
        ends = [p + 1 for p in positions]
        values = [d[p] * scale for p in positions]
        bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)

    bw.close()
    if log_path:
        log_fh.close()

if __name__ == "__main__":
    main()
