"""Synthetic checks for RNA orientation, unbiased selection, and motif statistics."""
import csv
import gzip
import importlib.util
import json
import io
from contextlib import redirect_stdout
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow/rules/scripts"
sys.path.insert(0, str(SCRIPTS))
from motif_utils import (bh_adjust, is_holdout, locus_key, locus_partition, motif_matches,
                         parse_meme_consensuses, revcomp, select_discovery_sequences)
from call_rt_termination_sites import load_counts, mapped_crosslinks, select_rt_sites, stream_selected_sites
from run_meme import run as run_meme
from motif_distribution import (add_pair, make_rows, matched_background, new_accumulator,
                                read_sites, run as run_distribution, sequence_at, load_blacklist, overlaps_blacklist)

try:
    import scipy
    import matplotlib
    import pysam
    SCIENTIFIC_RUNTIME = True
except ImportError:
    SCIENTIFIC_RUNTIME = False


class MotifUnitTests(unittest.TestCase):
    def test_rna_orientation_offset_and_multiple_matches(self):
        # Two DRACH motifs centered at 5 and 13 relative to a crosslink at 9.
        sequence = "CCCTGACTCCCTGACTCCC"
        self.assertEqual(motif_matches(sequence, 9, "DRACH", 3, 6), {-4, 4})
        # Reverse alignment is reverse-complemented back into the same RNA.
        self.assertEqual(motif_matches(revcomp(revcomp(sequence)), 9, "DRACH", 3, 6), {-4, 4})
        # All overlapping instances are represented; do not choose nearest one.
        self.assertEqual(motif_matches("AAAAA", 2, "AAA", 2, 1), {-1, 0, 1})
        self.assertEqual(motif_matches("CCCCCC", 3, "DRACH", 3, 2), set())
        with self.assertRaises(ValueError):
            motif_matches(sequence, 9, "DXACH", 3, 6)

    def test_site_selection_strict_threshold_neighbor_and_strand(self):
        counts = {("chr1", "+"): {100: 40, 200: 41, 300: 45, 299: 30, 400: 41, 401: 30},
                  ("chr1", "-"): {200: 1, 500: 50}}
        selected = select_rt_sites(counts, 40, 1.5)
        self.assertEqual(set(selected), {("chr1", 200, "+"), ("chr1", 300, "+"), ("chr1", 500, "-")})

    def test_streaming_call_matches_counts_across_strands_and_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "counts.bed"
            path.write_text("chr1\t1\t2\t.\t30\t+\nchr1\t1\t2\t.\t5\t-\n"
                            "chr1\t2\t3\t.\t45\t+\nchr1\t3\t4\t.\t10\t+\n"
                            "chr1\t9\t10\t.\t41\t-\nchr2\t0\t1\t.\t50\t-\n")
            expected = select_rt_sites(load_counts(path), 40, 1.5)
            observed, stats = stream_selected_sites(path, 40, 1.5)
            self.assertEqual(observed, expected)
            self.assertEqual(stats["distinct_r2_positions"], 6)
            self.assertEqual(stats["usable_read2"], 181)

    def test_mapped_splice_branches_and_unique_crosslinks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mapping.tsv.gz"
            with gzip.open(path, "wt") as out:
                out.write("chrom\tr2_pos0\tcrosslink_pos0\trna_strand\tcount\n"
                          "chr1\t200\t199\t+\t20\nchr1\t200\t99\t+\t21\n"
                          "chr1\t300\t301\t-\t41\n")
            selected = {("chr1", 200, "+"): (41, 0, 0), ("chr1", 300, "-"): (41, 0, 0)}
            sites, stats = mapped_crosslinks(path, selected)
            self.assertEqual(sites[("chr1", 99, "+")]["count"], 21)
            self.assertEqual(stats["splice_adjusted_reads"], 21)
            self.assertEqual(stats["crosslink_supporting_reads"], 82)
            self.assertEqual(stats["selected_r2_reads_without_valid_crosslink"], 0)

    def test_counts_are_unique_and_not_reads_as_replicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "count.bed"
            path.write_text("chr1\t5\t6\t.\t10\t+\n")
            self.assertEqual(load_counts(path)[("chr1", "+")][5], 10)
            self.assertEqual(read_sites(path)[("chr1", 5, "+")], 10)
            path.write_text(path.read_text() * 2)
            with self.assertRaises(ValueError):
                load_counts(path)
            with self.assertRaises(ValueError):
                read_sites(path)
        acc = new_accumulator(2)
        add_pair(acc, {-1, 1}, {1, 2}, 100, 2)
        self.assertEqual(acc["sites"], 1)
        self.assertEqual(acc["reads"], 100)
        self.assertEqual(acc["foreground"], [0, 1, 0, 1, 0])
        self.assertEqual(acc["foreground_only"], [0, 1, 0, 0, 0])
        self.assertEqual(acc["background_only"], [0, 0, 0, 0, 1])
        self.assertEqual(acc["foreground_multiple"], 1)

    def test_discovery_dedup_bounded_and_holdout_not_leaked(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "first.fa", Path(tmp) / "second.fa"
            first.write_text("".join(f">a{i}|locus=chr1:{i * 1000 + 500}:+\nACGTACGT\n" for i in range(200)))
            second.write_text(first.read_text())
            records, stats = select_discovery_sequences([first, second], 15, .2, 17)
            reverse, _ = select_discovery_sequences([second, first], 15, .2, 17)
            self.assertEqual(records, reverse)
            self.assertEqual(len(records), 15)
            self.assertEqual(stats["duplicate_loci"], 200)
            self.assertGreater(stats["holdout_loci"], 0)
            for name, seq in records:
                self.assertFalse(is_holdout(name.split("locus=")[1], .2, 17))
            # Same sequence at different sites is not collapsed.
            all_records, _ = select_discovery_sequences([first], 300, 0, 17)
            self.assertEqual(len(all_records), 200)

    def test_holdout_block_guards_prevent_neighbor_and_strand_leakage(self):
        self.assertEqual(locus_partition("chr1:500:+"), locus_partition("chr1:502:-"))
        self.assertEqual(locus_partition("chr1:9999:+"), "boundary_guard")
        self.assertEqual(locus_partition("chr1:10005:-"), "boundary_guard")
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "first.fa", Path(tmp) / "second.fa"
            first.write_text(">a|locus=chr1:500:+|center=r2_start:501\nACGTACGT\n"
                             ">b|locus=chr1:1500:+|center=r2_start:1501\nACGTACGT\n")
            second.write_text(">a2|locus=chr1:500:+|center=r2_start:600\nCCCCCCCC\n")
            one, stats = select_discovery_sequences([first, second], 10, 0, 17)
            two, _ = select_discovery_sequences([second, first], 10, 0, 17)
            self.assertEqual(one, two)
            self.assertEqual(len(one), 1)
            self.assertEqual(stats["ambiguous_duplicate_loci"], 1)

    def test_blacklist_interval_boundaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "blacklist.bed"
            path.write_text("chr1\t20\t30\nchr1\t25\t40\n")
            blacklist = load_blacklist([path])
            self.assertFalse(overlaps_blacklist(blacklist, "chr1", 10, 20))
            self.assertTrue(overlaps_blacklist(blacklist, "chr1", 10, 21))
            self.assertTrue(overlaps_blacklist(blacklist, "chr1", 39, 45))
            self.assertFalse(overlaps_blacklist(blacklist, "chr1", 40, 45))

    def test_bh_and_meme_consensus(self):
        adjusted = bh_adjust([.01, .04, .03, None])
        for observed, expected in zip(adjusted[:3], [.03, .04, .04]):
            self.assertAlmostEqual(observed, expected)
        self.assertIsNone(adjusted[3])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meme.txt"
            path.write_text("MOTIF 1 MEME-1\nletter-probability matrix: alength= 4 w= 3\n"
                            "0.9 0.04 0.03 0.03\n0.45 0.05 0.45 0.05\n0.05 0.9 0.03 0.02\n")
            model = parse_meme_consensuses(path)[0]
            self.assertEqual(model["sequence"], "ARC")
            self.assertEqual(model["anchor_pos"], 2)
            path.write_text("MEME not run: zero sites\n")
            self.assertEqual(parse_meme_consensuses(path), [])

    @unittest.skipUnless(SCIENTIFIC_RUNTIME, "pysam, scipy, matplotlib required")
    def test_paired_exact_test_and_equivalent_coordinates(self):
        acc = new_accumulator(1)
        for i in range(10):
            add_pair(acc, {0}, set(), 100, 1)
        motif = {"sequence": "DRACH", "anchor_pos": 3, "source": "predefined"}
        rows = make_rows("DRACH", motif, acc, 1, 1, r2_subset=True)
        center = rows[1]
        self.assertAlmostEqual(center["pvalue_one_sided_paired"], 1 / 1024)
        self.assertEqual(center["effective_site_pairs"], 10)
        self.assertEqual(center["supporting_reads"], 1000)
        self.assertEqual(center["crosslink_offset_from_motif_anchor_rna"], 0)
        self.assertEqual(center["r2_start_offset_from_motif_anchor_rna"], 1)
        empty = make_rows("DRACH", motif, new_accumulator(1), 1, 1)
        self.assertIsNone(empty[0]["pvalue_one_sided_paired"])
        self.assertIsNone(empty[0]["foreground_site_fraction"])

    @unittest.skipUnless(SCIENTIFIC_RUNTIME, "pysam, scipy, matplotlib required")
    def test_fasta_orientation_boundaries_and_repeatable_background(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reference.fa"
            seq = "C" * 500 + "CTGACT" + "C" * 1500
            path.write_text(">chr1\n" + seq + "\n")
            pysam.faidx(str(path))
            with pysam.FastaFile(str(path)) as fasta:
                lengths = dict(zip(fasta.references, fasta.lengths))
                forward, _ = sequence_at(fasta, "chr1", 503, "+", 4, lengths)
                reverse, _ = sequence_at(fasta, "chr1", 503, "-", 4, lengths)
                self.assertEqual(forward, revcomp(reverse))
                self.assertEqual(sequence_at(fasta, "chr1", 2, "+", 4, lengths), (None, "boundary"))
                bg1 = matched_background(fasta, "chr1", 503, "+", 4, lengths, 20, 100, 17, 20)
                bg2 = matched_background(fasta, "chr1", 503, "+", 4, lengths, 20, 100, 17, 20)
                self.assertEqual(bg1, bg2)
                self.assertGreaterEqual(abs(bg1[1] - 503), 20)

    @unittest.skipUnless(SCIENTIFIC_RUNTIME, "pysam, scipy, matplotlib required")
    def test_nonempty_distribution_strands_splice_subset_and_denominators(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, bed = root / "ref.fa", root / "sites.bed"
            sequence = list("C" * 10000)
            sequence[998:1003] = "TGACT"  # RNA + anchor at crosslink +0
            sequence[1996:2001] = "AGTCA"  # RNA - anchor at crosslink +2
            sequence[2997:3002] = "TGACT"  # RNA + anchor at crosslink -1
            sequence[4000] = "N"
            fasta.write_text(">chr1\n" + "".join(sequence) + "\n")
            pysam.faidx(str(fasta))
            bed.write_text("chr1\t1000\t1001\ta\t10\t+\t1001\n"
                           "chr1\t2000\t2001\tb\t20\t-\t1999\n"
                           "chr1\t3000\t3001\tc\t30\t+\t3100\n"
                           "chr1\t2\t3\td\t5\t+\t3\n"
                           "chr1\t4000\t4001\te\t5\t+\t4001\n")
            smk = SimpleNamespace(wildcards=SimpleNamespace(sample="synthetic"),
                  params=SimpleNamespace(window=5, seed=17, holdout_fraction=.2, site_shift=1,
                         holdout_block_size=10000, holdout_guard=100,
                         motifs=[{"name": "DRACH", "sequence": "DRACH", "anchor_pos": 3}],
                         consensus_mass=.8, background_min_distance=200, background_max_distance=1000,
                         background_attempts=20),
                  input=SimpleNamespace(bed=str(bed), fasta=str(fasta), meme=[]),
                  output=SimpleNamespace(table=str(root / "offsets.tsv"), r2_table=str(root / "r2_offsets.tsv"), stats=str(root / "stats.json"),
                         png=str(root / "plot.png"), pdf=str(root / "plot.pdf"), pairs=str(root / "pairs.tsv")))
            with redirect_stdout(io.StringIO()):
                run_distribution(smk)
            report = json.loads((root / "stats.json").read_text())
            self.assertEqual(report["motifs"][0]["effective_site_pairs"], 3)
            self.assertEqual(report["r2_contiguous_subset_motifs"][0]["effective_site_pairs"], 2)
            self.assertEqual(report["counts"]["foreground_boundary_skipped"], 1)
            self.assertEqual(report["counts"]["foreground_ambiguous_bases_skipped"], 1)
            with (root / "offsets.tsv").open() as handle:
                rows = {int(row["motif_anchor_offset_from_crosslink_rna"]): row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual([int(rows[d]["foreground_sites"]) for d in [-1, 0, 2]], [1, 1, 1])
            self.assertEqual(rows[0]["r2_start_offset_from_motif_anchor_rna"], "NA")
            with (root / "r2_offsets.tsv").open() as handle:
                rows = {int(row["r2_start_offset_from_motif_anchor_rna"]): row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(int(rows[1]["foreground_sites"]), 1)
            self.assertEqual(int(rows[-1]["foreground_sites"]), 1)
            self.assertEqual(int(rows[2]["foreground_sites"]), 0)

    @unittest.skipUnless(SCIENTIFIC_RUNTIME, "pysam, scipy, matplotlib required")
    def test_empty_meme_emits_decodable_png_and_explicit_stats(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "empty.fa"
            fasta.write_text("")
            outputs = {key: str(root / filename) for key, filename in
                       (("html", "meme.html"), ("txt", "meme.txt"), ("logo", "logo.png"),
                        ("sequences", "discovery.fa"), ("stats", "stats.json"))}
            smk = SimpleNamespace(input=SimpleNamespace(fasta=str(fasta)), output=outputs,
                  params=SimpleNamespace(label="empty", min_sequences=10, max_sequences=100,
                         holdout_fraction=.2, seed=17, holdout_block_size=10000, holdout_guard=100),
                  log=[str(root / "meme.log")])
            run_meme(smk)
            with Image.open(outputs["logo"]) as logo:
                logo.load()
                self.assertEqual(logo.size, (1, 1))
            self.assertEqual(json.loads(Path(outputs["stats"]).read_text())["status"], "insufficient_discovery_sequences")

    @unittest.skipUnless(SCIENTIFIC_RUNTIME, "pysam, scipy, matplotlib required")
    def test_empty_distribution_complete_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, bed = root / "ref.fa", root / "sites.bed"
            fasta.write_text(">chr1\n" + "ACGT" * 100 + "\n")
            pysam.faidx(str(fasta))
            bed.write_text("")
            smk = SimpleNamespace(wildcards=SimpleNamespace(sample="synthetic_empty"),
                  params=SimpleNamespace(window=5, seed=17, holdout_fraction=.2, site_shift=1, holdout_block_size=10000, holdout_guard=100,
                         motifs=[{"name": "DRACH", "sequence": "DRACH", "anchor_pos": 3}],
                         consensus_mass=.8, background_min_distance=30, background_max_distance=100,
                         background_attempts=5),
                  input=SimpleNamespace(bed=str(bed), fasta=str(fasta), meme=[]),
                  output=SimpleNamespace(table=str(root / "offsets.tsv"), r2_table=str(root / "r2_offsets.tsv"), stats=str(root / "stats.json"),
                         png=str(root / "plot.png"), pdf=str(root / "plot.pdf"), pairs=str(root / "pairs.tsv")))
            with redirect_stdout(io.StringIO()):
                run_distribution(smk)
            report = json.loads((root / "stats.json").read_text())
            self.assertEqual(report["motifs"][0]["status"], "no_effective_sites")
            self.assertEqual(report["motifs"][0]["peak_frequency_offsets_nt"], [])
            with (root / "offsets.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 11)
            self.assertEqual(rows[0]["qvalue_bh_sample"], "NA")
            self.assertGreater((root / "plot.png").stat().st_size, 100)
            self.assertGreater((root / "plot.pdf").stat().st_size, 100)


if __name__ == "__main__":
    unittest.main()
