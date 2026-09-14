"""Count unbiased Read-2 endpoints and CIGAR-aware inferred cross-link sites.

BED6 strand is the RNA strand (opposite Read 2 in this GRIP library).
A cross-link is inferred one nucleotide upstream of an endpoint in RNA space;
this is not a motif-conditioned or confirmed modification coordinate. Five-prime
unaligned sequence is excluded by default because its RT boundary is unknown.
Coordinate-sorted BAMs are counted with a moving heap, not a genome-sized dict.
"""
import argparse
import gzip
import heapq
import json
import os
import sys
import tempfile
from bisect import bisect_right
from collections import Counter
from contextlib import ExitStack
from types import SimpleNamespace

import pysam
import pyBigWig

MATCH = {0, 7, 8}
REF = {0, 2, 3, 7, 8}


def load_fai(path):
    with open(path) as handle:
        chroms = [(f[0], int(f[1])) for line in handle if (f := line.split())]
    if len(dict(chroms)) != len(chroms) or any(n <= 0 for _, n in chroms):
        raise ValueError("FAI must have unique chromosomes with positive lengths")
    return chroms


def load_blacklist(path):
    """Merge overlapping/nested BED intervals before binary-search membership."""
    intervals = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError(f"Invalid blacklist interval: {line.rstrip()}")
            intervals.setdefault(chrom, []).append((start, end))
    merged = {}
    for chrom, entries in intervals.items():
        result = []
        for start, end in sorted(entries):
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        merged[chrom] = ([s for s, _ in result], result)

    def contains(chrom, pos):
        if chrom not in merged:
            return False
        starts, entries = merged[chrom]
        i = bisect_right(starts, pos) - 1
        return i >= 0 and pos < entries[i][1]

    return contains


def endpoint_and_crosslink(read, shift=1):
    """Return first aligned query base, inferred RNA-upstream site and flags.

    Walk from the sequencing 5-prime end. Terminal D/N are not query bases.
    For this opposite-strand library RNA-upstream advances into the R2 CIGAR.
    N is an intron and consumes no transcript bases; D consumes a reference
    transcript base. Soft/hard clips and insertions before the first aligned
    base are reported explicitly, never silently extrapolated onto the genome.
    """
    if shift < 0:
        raise ValueError("crosslink_shift must be nonnegative")
    cigar = read.cigartuples or []
    direction = -1 if read.is_reverse else 1
    cursor = read.reference_end - 1 if read.is_reverse and read.reference_end is not None else read.reference_start
    if read.is_reverse:
        cigar = reversed(cigar)
    flags = set()
    endpoint = None
    remaining = shift
    for op, length in cigar:
        if endpoint is None:
            if op in MATCH:
                endpoint = cursor
            elif op in (4, 5, 1):
                flags.add({4: "softclip", 5: "hardclip", 1: "insertion"}[op])
            elif op in (2, 3):
                flags.add("terminal_reference_gap")
        if endpoint is not None and op in MATCH | {2}:
            if remaining < length:
                return endpoint, cursor + direction * remaining, flags
            remaining -= length
        if op in REF:
            cursor += direction * length
    return endpoint, None, flags


class SlidingCounts:
    """Keep only coordinates ahead of the sorted BAM reference-start watermark."""
    def __init__(self, emit):
        self.heap = []
        self.counts = {}
        self.emit = emit
        self.max_pending = 0

    def add(self, key):
        if key not in self.counts:
            heapq.heappush(self.heap, key)
            self.counts[key] = 0
        self.counts[key] += 1
        self.max_pending = max(self.max_pending, len(self.counts))

    def flush(self, before=float("inf")):
        while self.heap and self.heap[0][0] < before:
            key = heapq.heappop(self.heap)
            self.emit(key, self.counts.pop(key))


def _open_text(path, mode):
    return gzip.open(path, mode + "t") if str(path).endswith(".gz") else open(path, mode)


def write_bigwigs(bed_path, paths, chroms, accepted):
    """Second streaming pass after the shared two-strand CPM denominator is known."""
    paths = {k: p for k, p in paths.items() if p}
    if not paths:
        return
    buffers = {k: ([], [], []) for k in paths}
    scale = 1e6 / accepted if accepted else 0.0
    with ExitStack() as stack:
        writers = {k: stack.enter_context(pyBigWig.open(str(p), "w")) for k, p in paths.items()}
        for writer in writers.values():
            writer.addHeader(chroms)

        def flush(key):
            names, starts, vals = buffers[key]
            if starts:
                writers[key].addEntries(names, starts, ends=[x + 1 for x in starts], values=vals)
                names.clear(); starts.clear(); vals.clear()

        def append(key, chrom, pos, count):
            if key in writers:
                names, starts, vals = buffers[key]
                names.append(chrom); starts.append(pos); vals.append(float(count * scale))
                if len(starts) >= 65536:
                    flush(key)

        pending = None
        total = 0
        with _open_text(bed_path, "r") as handle:
            for line in handle:
                chrom, start, _, _, count, strand = line.split()[:6]
                pos, count = int(start), int(count)
                if pending != (chrom, pos):
                    if pending is not None:
                        append("bw", *pending, total)
                    pending, total = (chrom, pos), 0
                total += count
                append("plus" if strand == "+" else "minus", chrom, pos, count)
            if pending is not None:
                append("bw", *pending, total)
        for key in writers:
            flush(key)


def extract(args):
    if args.five_prime_clip not in ("exclude", "aligned"):
        raise ValueError("five_prime_clip must be exclude or aligned")
    if args.crosslink_shift < 0 or not 0 <= args.min_mapq <= 255:
        raise ValueError("crosslink_shift must be nonnegative and min_mapq within 0..255")
    chroms_fai = load_fai(args.fai)
    lengths = dict(chroms_fai)
    excluded = load_blacklist(args.blacklist_bed) if args.blacklist_bed else lambda c, p: False
    stats = Counter()
    sample = getattr(args, "sample", "site")
    with tempfile.TemporaryDirectory(prefix="grip-signal-") as tmp, ExitStack() as stack:
        endpoint_bed = args.out_bed or os.path.join(tmp, "endpoints.bed.gz")
        crosslink_bed = getattr(args, "out_crosslink_bed", None) or os.path.join(tmp, "crosslinks.bed.gz")
        mapping_tsv = getattr(args, "out_mapping", None) or os.path.join(tmp, "mapping.tsv.gz")
        paths = [endpoint_bed, crosslink_bed, mapping_tsv]
        for name in ("out_bw", "out_plus", "out_minus", "out_crosslink_bw", "out_crosslink_plus", "out_crosslink_minus", "out_stats"):
            if getattr(args, name, None):
                paths.append(getattr(args, name))
        for path in paths:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        bed = stack.enter_context(_open_text(endpoint_bed, "w"))
        crossbed = stack.enter_context(_open_text(crosslink_bed, "w"))
        mapping = stack.enter_context(_open_text(mapping_tsv, "w"))
        mapping.write("chrom\tr2_pos0\tcrosslink_pos0\trna_strand\tcount\n")
        bam = stack.enter_context(pysam.AlignmentFile(args.bam, "rb"))
        chroms = [(c, lengths[c]) for c in bam.references if c in lengths]
        chroms.extend((c, n) for c, n in chroms_fai if c not in bam.references)
        for chrom, length in zip(bam.references, bam.lengths):
            if chrom in lengths and length != lengths[chrom]:
                raise ValueError(f"BAM/FAI chromosome length mismatch for {chrom}")
        chrom = None

        def emit_bed(handle, label, key, count):
            pos, strand = key
            handle.write(f"{chrom}\t{pos}\t{pos + 1}\t{sample}_{label}_{chrom}_{pos + 1}_{strand}\t{count}\t{strand}\n")
            stats[label + "_unique_sites"] += 1

        endpoints = SlidingCounts(lambda key, n: emit_bed(bed, "r2", key, n))
        crosslinks = SlidingCounts(lambda key, n: emit_bed(crossbed, "crosslink", key, n))
        pairs = SlidingCounts(lambda key, n: mapping.write(f"{chrom}\t{key[0]}\t{key[1]}\t{key[2]}\t{n}\n"))
        counters = (endpoints, crosslinks, pairs)
        previous = (-1, -1)
        for read in bam.fetch(until_eof=True):
            stats["records_total"] += 1
            if read.is_unmapped:
                stats["unmapped"] += 1
                continue
            coordinate = (read.reference_id, read.reference_start)
            if coordinate < previous:
                raise ValueError("Input BAM must be coordinate sorted (out-of-order record found)")
            previous = coordinate
            next_chrom = bam.get_reference_name(read.reference_id)
            if next_chrom != chrom:
                for counter in counters:
                    counter.flush()
                chrom = next_chrom
            else:
                for counter in counters:
                    counter.flush(read.reference_start)
            if not read.is_read2 or not read.is_paired:
                continue
            stats["read2_total"] += 1
            reason = None
            if read.is_secondary or read.is_supplementary:
                reason = "nonprimary"
            elif getattr(args, "exclude_qcfail", True) and read.is_qcfail:
                reason = "qcfail"
            elif getattr(args, "exclude_duplicates", False) and read.is_duplicate:
                reason = "duplicate"
            elif read.mapping_quality < args.min_mapq:
                reason = "low_mapq"
            elif args.require_proper_pair and not read.is_proper_pair:
                reason = "improper_pair"
            elif chrom not in lengths:
                reason = "unknown_reference"
            if reason:
                stats["excluded_" + reason] += 1
                continue
            stats["read2_pre_endpoint_filter"] += 1
            endpoint, crosslink, flags = endpoint_and_crosslink(read, args.crosslink_shift)
            for flag in flags:
                stats["five_prime_" + flag] += 1
            unaligned = bool(flags & {"softclip", "hardclip", "insertion"})
            if unaligned:
                # A read can have H, S and I together; the QC fraction needs a union.
                stats["five_prime_unaligned_any"] += 1
            if args.five_prime_clip == "exclude" and unaligned:
                stats["excluded_five_prime_unaligned"] += 1
                continue
            if endpoint is None or crosslink is None:
                stats["excluded_no_aligned_endpoint_or_shift"] += 1
                continue
            if not (0 <= endpoint < lengths[chrom] and 0 <= crosslink < lengths[chrom]):
                stats["excluded_reference_boundary"] += 1
                continue
            if excluded(chrom, endpoint) or excluded(chrom, crosslink):
                stats["excluded_blacklist_endpoint_or_crosslink"] += 1
                continue
            strand = "+" if read.is_reverse else "-"
            endpoints.add((endpoint, strand))
            crosslinks.add((crosslink, strand))
            pairs.add((endpoint, crosslink, strand))
            stats["accepted_read2"] += 1
            stats["accepted_rna_" + ("plus" if strand == "+" else "minus")] += 1
            if abs(crosslink - endpoint) != args.crosslink_shift:
                stats["splice_aware_shift"] += 1
        for counter in counters:
            counter.flush()
        for handle in (bed, crossbed, mapping):
            handle.close()
        stats["maximum_pending_coordinates"] = sum(c.max_pending for c in counters)
        n = stats["accepted_read2"]
        write_bigwigs(endpoint_bed, {"bw": args.out_bw, "plus": getattr(args, "out_plus", None), "minus": getattr(args, "out_minus", None)}, chroms, n)
        write_bigwigs(crosslink_bed, {"bw": getattr(args, "out_crosslink_bw", None), "plus": getattr(args, "out_crosslink_plus", None), "minus": getattr(args, "out_crosslink_minus", None)}, chroms, n)
    result = dict(sorted(stats.items()))
    result.update({"cpm_denominator": stats["accepted_read2"], "cpm_denominator_definition": "accepted primary mapped Read2 records after endpoint, clipping, pair, QC, MAPQ and joint endpoint/crosslink blacklist filtering; shared by both RNA strands", "rna_strand": "opposite_read2", "coordinate_system": "0-based half-open BED; 0-based TSV positions", "crosslink_shift_rna_upstream": args.crosslink_shift, "five_prime_clip_policy": args.five_prime_clip, "empty_signal": stats["accepted_read2"] == 0})
    if getattr(args, "out_stats", None):
        with open(args.out_stats, "w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    print(json.dumps(result, indent=2))
    return result


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fai", required=True)
    for name in ("out-bw", "out-bed", "out-plus", "out-minus", "out-crosslink-bw", "out-crosslink-bed", "out-crosslink-plus", "out-crosslink-minus", "out-mapping", "out-stats", "blacklist-bed"):
        parser.add_argument("--" + name)
    parser.add_argument("--sample", default="site")
    parser.add_argument("--min-mapq", type=int, default=11)
    parser.add_argument("--crosslink-shift", type=int, default=1)
    parser.add_argument("--five-prime-clip", choices=("exclude", "aligned"), default="exclude")
    parser.add_argument("--require-proper-pair", action="store_true")
    parser.add_argument("--exclude-duplicates", action="store_true")
    parser.add_argument("--include-qcfail", dest="exclude_qcfail", action="store_false")
    return parser.parse_args()


def resolve_args():
    if "snakemake" not in globals():
        return parse_args(), None
    output_names = {"out_bw": "bw", "out_bed": "bed", "out_plus": "plus", "out_minus": "minus", "out_crosslink_bw": "crosslink_bw", "out_crosslink_bed": "crosslink_bed", "out_crosslink_plus": "crosslink_plus", "out_crosslink_minus": "crosslink_minus", "out_mapping": "mapping", "out_stats": "stats"}
    args = {key: snakemake.output.get(value) for key, value in output_names.items()}
    args.update(bam=snakemake.input.bam, fai=snakemake.input.fai, blacklist_bed=(str(snakemake.input.bl[0]) if snakemake.input.get("bl") else None), sample=str(snakemake.wildcards.sample), min_mapq=int(snakemake.params.min_mapq), crosslink_shift=int(snakemake.params.get("crosslink_shift", 1)), five_prime_clip=str(snakemake.params.get("five_prime_clip", "exclude")), require_proper_pair=bool(snakemake.params.get("require_proper_pair", False)), exclude_duplicates=bool(snakemake.params.get("exclude_duplicates", False)), exclude_qcfail=bool(snakemake.params.get("exclude_qcfail", True)))
    return SimpleNamespace(**args), str(snakemake.log[0]) if snakemake.log else None


if __name__ == "__main__":
    args, log_path = resolve_args()
    if log_path:
        os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
        with open(log_path, "w") as log:
            from contextlib import redirect_stdout, redirect_stderr
            with redirect_stdout(log), redirect_stderr(log):
                extract(args)
    else:
        extract(args)
