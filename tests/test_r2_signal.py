"""Synthetic read-2 tests exercise coordinate, strand and denominator invariants."""
import contextlib
import gzip
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pysam
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2_signal", ROOT / "workflow/rules/scripts/extract_r2_first_base.py")
signal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(signal)


def read(start, cigar, reverse=False, extra_flags=0, mapq=60, name="r"):
    r = pysam.AlignedSegment()
    r.query_name = name
    r.flag = 1 | 2 | 128 | (16 if reverse else 0) | extra_flags
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = mapq
    r.cigarstring = cigar
    r.query_sequence = "A" * sum(n for op, n in r.cigartuples if op in (0, 1, 4, 7, 8))
    return r


class EndpointTests(unittest.TestCase):
    def test_both_alignment_orientations(self):
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "20M"))[:2], (10, 9))
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "20M", True))[:2], (29, 30))

    def test_terminal_reference_gaps_not_query_bases(self):
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "3D20M"))[:2], (13, 12))
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "20M3D", True))[:2], (29, 30))
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "3N20M"))[:2], (13, 12))
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "20M3N", True))[:2], (29, 30))

    def test_crosslink_shift_skips_introns_on_both_strands(self):
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "1M100N19M"), rna_strand="opposite_read2")[:2], (10, 111))
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "19M100N1M", True), rna_strand="opposite_read2")[:2], (129, 28))

    def test_deletion_retains_reference_transcript_bases(self):
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "1M2D19M"))[:2], (10, 9))
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "19M2D1M", True))[:2], (31, 32))

    def test_five_prime_clipping_and_insertion_flags(self):
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "3S20M"))[2], {"softclip"})
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "20M3H", True))[2], {"hardclip"})
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "2I20M"))[2], {"insertion"})
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "20M3S"))[2], set())
        self.assertEqual(signal.endpoint_and_crosslink(read(10, "3H20M", True))[2], set())

    def test_short_alignment_not_extrapolated_past_boundary(self):
        self.assertEqual(signal.endpoint_and_crosslink(read(0, "1M"))[:2], (0, -1))
        self.assertEqual(signal.endpoint_and_crosslink(read(0, "1M"), 0)[:2], (0, 0))
        with self.assertRaises(ValueError):
            signal.endpoint_and_crosslink(read(0, "1M"), -1)

    def test_nested_blacklist_and_half_open_boundaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "blacklist.bed"
            bed.write_text("chr1\t10\t100\nchr1\t20\t30\nchr1\t95\t120\n")
            contains = signal.load_blacklist(bed)
            self.assertTrue(contains("chr1", 50))
            self.assertTrue(contains("chr1", 119))
            self.assertFalse(contains("chr1", 120))
            self.assertFalse(contains("chr1", 9))


class ExtractionTests(unittest.TestCase):
    def run_signal(self, tmp, reads, **options):
        tmp = Path(tmp)
        bam = tmp / "input.bam"
        with pysam.AlignmentFile(bam, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}) as handle:
            for r in reads:
                handle.write(r)
        fai = tmp / "reference.fa.fai"
        fai.write_text("chr1\t1000\t0\t0\t0\n")
        args = dict(bam=str(bam), fai=str(fai), out_bw=str(tmp / "r2.bw"), out_plus=str(tmp / "r2.plus.bw"), out_minus=str(tmp / "r2.minus.bw"), out_crosslink_bw=str(tmp / "cross.bw"), out_crosslink_plus=str(tmp / "cross.plus.bw"), out_crosslink_minus=str(tmp / "cross.minus.bw"), out_bed=str(tmp / "r2.bed.gz"), out_crosslink_bed=str(tmp / "cross.bed.gz"), out_mapping=str(tmp / "mapping.tsv.gz"), out_stats=str(tmp / "stats.json"), min_mapq=11, blacklist_bed=None, require_proper_pair=False, crosslink_shift=1, five_prime_clip="exclude", exclude_qcfail=True, exclude_duplicates=False, sample="synthetic")
        args.update(options)
        with contextlib.redirect_stdout(io.StringIO()):
            result = signal.extract(SimpleNamespace(**args))
        return result, args

    def test_filters_cpm_and_rna_strand_contract(self):
        reads = [read(10, "20M"), read(10, "20M"), read(10, "20M", True), read(10, "2S20M"), read(10, "20M", extra_flags=512), read(10, "20M", extra_flags=256), read(10, "20M", extra_flags=2048), read(10, "20M", mapq=10), read(10, "20M", extra_flags=1024)]
        with tempfile.TemporaryDirectory() as tmp:
            stats, args = self.run_signal(tmp, reads)
            self.assertEqual(stats["accepted_read2"], 4)  # Dup flag opt-in, not UMI dedup.
            self.assertEqual(stats["cpm_denominator"], 4)
            self.assertEqual(stats["excluded_five_prime_unaligned"], 1)
            self.assertEqual(stats["excluded_nonprimary"], 2)
            self.assertEqual(stats["excluded_qcfail"], 1)
            with pyBigWig.open(args["out_bw"]) as bw:
                self.assertEqual(bw.values("chr1", 10, 11)[0], 750000)
                self.assertEqual(bw.values("chr1", 29, 30)[0], 250000)
            with pyBigWig.open(args["out_plus"]) as bw:
                self.assertEqual(bw.values("chr1", 10, 11)[0], 750000)
            with pyBigWig.open(args["out_minus"]) as bw:
                self.assertEqual(bw.values("chr1", 29, 30)[0], 250000)
            with gzip.open(args["out_bed"], "rt") as f:
                rows = [l.split() for l in f]
            self.assertEqual([(int(r[1]), int(r[4]), r[5]) for r in rows], [(10, 3, "+"), (29, 1, "-")])

    def test_splice_branch_mapping_counts_conserve_endpoint_total(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats, args = self.run_signal(tmp, [read(10, "20M"), read(10, "1M100N19M"), read(10, "1M100N19M"), read(30, "20M")], rna_strand="opposite_read2")
            self.assertEqual(stats["splice_aware_shift"], 2)
            with gzip.open(args["out_mapping"], "rt") as f:
                rows = [l.strip().split("\t") for l in f][1:]
            self.assertEqual(rows, [["chr1", "10", "11", "-", "1"], ["chr1", "10", "111", "-", "2"], ["chr1", "30", "31", "-", "1"]])
            with pyBigWig.open(args["out_crosslink_bw"]) as bw:
                self.assertEqual(bw.values("chr1", 111, 112)[0], 500000)

    def test_joint_endpoint_crosslink_blacklist_and_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            bl = Path(tmp) / "bl.bed"
            bl.write_text("chr1\t9\t10\n")  # Only inferred crosslink overlaps.
            stats, args = self.run_signal(tmp, [read(10, "20M"), read(20, "20M")], blacklist_bed=str(bl))
            self.assertEqual(stats["cpm_denominator"], 1)
            self.assertEqual(stats["excluded_blacklist_endpoint_or_crosslink"], 1)
            with pyBigWig.open(args["out_bw"]) as bw:
                self.assertEqual(bw.values("chr1", 20, 21)[0], 1000000)

    def test_upstream_shift_preserves_sorted_crosslink_bed(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats, args = self.run_signal(tmp, [read(10, "1M", True), read(11, "1M", True), read(12, "1M"), read(12, "2M")])
            with gzip.open(args["out_crosslink_bed"], "rt") as handle:
                rows = [line.split() for line in handle]
            self.assertEqual([(int(r[1]), r[5], int(r[4])) for r in rows], [(11, "+", 2), (11, "-", 1), (12, "-", 1)])
            self.assertEqual(stats["rna_strand"], "same_as_read2")

    def test_reference_boundaries_exclude_upstream_crosslinks(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats, _ = self.run_signal(tmp, [read(0, "1M"), read(10, "1M"), read(999, "1M", True)])
            self.assertEqual(stats["accepted_read2"], 1)
            self.assertEqual(stats["excluded_reference_boundary"], 2)

    def test_published_library_motif_coordinates_on_both_strands(self):
        # Independent fixture: known RNA GGACT has its central A at +3 from
        # the physical crosslink; R2 begins one base downstream of crosslink.
        genome = list("C" * 1000)
        genome[11:16] = "GGACT"  # RNA A at 13; crosslink 10, forward R2 11.
        genome[184:189] = "AGTCC"  # RNA A at 186; crosslink 189, reverse R2 188.
        with tempfile.TemporaryDirectory() as tmp:
            _, args = self.run_signal(tmp, [read(11, "20M"), read(169, "20M", True)])
            with gzip.open(args["out_crosslink_bed"], "rt") as handle:
                rows = [line.split() for line in handle]
            self.assertEqual([(int(r[1]), r[5]) for r in rows], [(10, "+"), (189, "-")])
            for row in rows:
                p, strand = int(row[1]), row[5]
                sequence = "".join(genome[p-5:p+6])
                if strand == "-":
                    sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                self.assertEqual(sequence[6:11], "GGACT")

    def test_aligned_sensitivity_policy_keeps_clipped_reads_with_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats, _ = self.run_signal(tmp, [read(10, "5S20M")], five_prime_clip="aligned")
            self.assertEqual(stats["five_prime_softclip"], 1)
            self.assertEqual(stats["accepted_read2"], 1)

    def test_mixed_five_prime_unaligned_operations_count_one_read(self):
        for policy in ("exclude", "aligned"):
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as tmp:
                stats, _ = self.run_signal(tmp, [read(10, "5H3S2I20M")], five_prime_clip=policy)
                self.assertEqual(stats["five_prime_softclip"], 1)
                self.assertEqual(stats["five_prime_hardclip"], 1)
                self.assertEqual(stats["five_prime_insertion"], 1)
                self.assertEqual(stats["five_prime_unaligned_any"], 1)
                self.assertEqual(stats["read2_pre_endpoint_filter"], 1)
                self.assertEqual(stats.get("accepted_read2", 0), int(policy == "aligned"))

    def test_empty_signal_is_explicit_and_writable(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats, args = self.run_signal(tmp, [read(10, "5S20M")])
            self.assertTrue(stats["empty_signal"])
            self.assertEqual(stats["cpm_denominator"], 0)
            with pyBigWig.open(args["out_bw"]) as bw:
                self.assertEqual(bw.chroms(), {"chr1": 1000})

    def test_unsorted_input_fails_instead_of_corrupting_tracks(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "coordinate sorted"):
                self.run_signal(tmp, [read(20, "20M"), read(10, "20M")])


class MotifDiagnosticTests(unittest.TestCase):
    def test_rna_strand_and_no_second_crosslink_shift(self):
        spec = importlib.util.spec_from_file_location("motif_anchor", ROOT / "workflow/rules/scripts/motif_anchor.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            seq = list("A" * 60)
            seq[10:15] = "GGACT"
            seq[30:35] = "AGTCC"  # Reverse complement of GGACT.
            fasta = tmp / "ref.fa"
            fasta.write_text(">chr1\n" + "".join(seq) + "\n")
            pysam.faidx(str(fasta))
            bed = tmp / "sites.bed"
            bed.write_text("chr1\t12\t13\tplus\t2\t+\nchr1\t32\t33\tminus\t3\t-\n")
            out = tmp / "annotated.bed"
            args = SimpleNamespace(bed_in=str(bed), fasta=str(fasta), out_bed=str(out), motif="DRACH", motif_pos=3, site_shift=0, search_window=0, keep_unmatched=False)
            with contextlib.redirect_stdout(io.StringIO()):
                result = module.annotate(args)
            self.assertEqual(result["matched_reads"], 5)
            rows = [x.split() for x in out.read_text().splitlines()]
            self.assertEqual([(r[1], r[5], r[6], r[8]) for r in rows], [("12", "+", "GGACT", "0"), ("32", "-", "GGACT", "0")])

    def test_equidistant_motif_hits_are_reported_ambiguous(self):
        spec = importlib.util.spec_from_file_location("motif_anchor", ROOT / "workflow/rules/scripts/motif_anchor.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = tmp / "ref.fa"
            fasta.write_text(">chr1\nAACACAA\n")
            pysam.faidx(str(fasta))
            bed = tmp / "sites.bed"
            bed.write_text("chr1\t3\t4\tplus\t2\t+\n")
            out = tmp / "annotated.bed"
            args = SimpleNamespace(bed_in=str(bed), fasta=str(fasta), out_bed=str(out), motif="C", motif_pos=1, site_shift=0, search_window=1, keep_unmatched=True)
            with contextlib.redirect_stdout(io.StringIO()):
                result = module.annotate(args)
            self.assertEqual(result["ambiguous_reads"], 2)
            self.assertEqual(out.read_text().split()[1], "3")
            self.assertEqual(out.read_text().split()[-1], "ambiguous")


if __name__ == "__main__":
    unittest.main()
