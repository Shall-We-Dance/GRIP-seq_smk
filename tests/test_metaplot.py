"""Analytic tests for mature-RNA geometry, sparse integration and density output."""

import csv
import json
from pathlib import Path
import runpy
import sys
import tempfile
from types import SimpleNamespace
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "rules" / "scripts"
sys.path.insert(0, str(SCRIPTS))
from metagene import common_chromosome_sizes, prepare_annotation, segment_transcript, scaled_region_profile

try:
    import numpy as np
except ImportError:
    np = None


class Named(dict):
    """Snakemake-like named IO lists for executing scripts without a scheduler."""
    def __getattr__(self, name):
        return self[name]

    def __iter__(self):
        return iter(self.values())


class FakeBigWig:
    def __init__(self, intervals):
        self.data = intervals

    def intervals(self, chrom, start, end):
        return [row for row in self.data if row[0] < end and row[1] > start]


def record(strand="+"):
    return {"gene_id": "g", "transcript_id": "t", "chrom": "chrTest", "strand": strand,
            "exons": [[100, 110], [200, 212]], "cds": [[105, 110], [200, 206]],
            "stop_codons": [[206, 209]]}


class AnnotationTests(unittest.TestCase):
    def test_shared_contigs_require_consistent_lengths(self):
        sizes, summary = common_chromosome_sizes([
            {"chrTest": 100, "chrPatch": 80}, {"chrTest": 100, "chrOther": 90},
        ])
        self.assertEqual(sizes, {"chrTest": 100})
        self.assertEqual(summary["contigs_not_shared"], ["chrOther", "chrPatch"])
        self.assertEqual(summary["track_headers_checked"], 2)
        with self.assertRaisesRegex(ValueError, "disagree on contig lengths"):
            common_chromosome_sizes([{"chrTest": 100}, {"chrTest": 101}])
        with self.assertRaisesRegex(ValueError, "share no common contigs"):
            common_chromosome_sizes([{"chrA": 100}, {"chrB": 100}])

    def test_patch_candidate_filtered_before_representative_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.gtf"
            lines = []
            for chrom, gene, tx, start, end, cds_start, cds_end in [
                ("chrPatch", "rescued", "preferred", 1, 40, 6, 30),
                ("chrTest", "rescued", "primary_2", 1, 30, 6, 20),
                ("chrPatch", "patch_only", "missing", 51, 80, 56, 70),
                ("chrTest", "out_of_bounds", "beyond_end", 91, 120, 96, 110),
            ]:
                for feature, a, b in [("exon", start, end), ("CDS", cds_start, cds_end)]:
                    lines.append(f'{chrom}\ttest\t{feature}\t{a}\t{b}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{tx}";\n')
            path.write_text("".join(lines))
            models, rows, summary = prepare_annotation(path, available_chrom_sizes={"chrTest": 100})
            self.assertEqual(len(models), 1)
            self.assertEqual(models[0]["transcript_id"], "primary_2")
            self.assertEqual(summary["genes_with_available_contig_fallback"], 1)
            self.assertEqual(summary["genes_no_compatible_contig_candidate"], 2)
            self.assertEqual(summary["transcripts_structurally_eligible"], 4)
            self.assertEqual(summary["transcripts_eligible"], 1)
            self.assertEqual(summary["transcript_exclusions"]["contig_not_shared_by_all_tracks"], 2)
            self.assertEqual(summary["transcript_exclusions"]["exons_outside_shared_contig_length"], 1)
            row = next(row for row in rows if row["gene_id"] == "rescued")
            self.assertEqual(row["reason"], "available_contig_fallback")
            self.assertEqual(row["preferred_chrom_before_contig_filter"], "chrPatch")
            self.assertEqual(summary["region_median_lengths"], {"utr5": 5, "cds": 15, "utr3": 10})

    def test_plus_exons_and_exclusive_stop_codon(self):
        model = segment_transcript(record())
        self.assertEqual(model["lengths"], {"utr5": 5, "cds": 14, "utr3": 3})
        self.assertEqual(model["regions"]["cds"], [[105, 110], [200, 209]])
        self.assertEqual(model["regions"]["utr3"], [[209, 212]])

    def test_minus_orientation_and_stop_codon(self):
        tx = record("-")
        tx.update(exons=[[100, 112], [200, 210]], cds=[[106, 112], [200, 205]], stop_codons=[[103, 106]])
        model = segment_transcript(tx)
        self.assertEqual(model["lengths"], {"utr5": 5, "cds": 14, "utr3": 3})
        self.assertEqual(model["regions"]["utr5"], [[205, 210]])
        self.assertEqual(model["regions"]["cds"], [[200, 205], [103, 112]])
        self.assertEqual(model["regions"]["utr3"], [[100, 103]])

    def test_stop_codon_split_over_junction(self):
        tx = record()
        tx.update(cds=[[105, 108]], stop_codons=[[108, 110], [200, 201]])
        model = segment_transcript(tx)
        self.assertEqual(model["regions"]["cds"], [[105, 110], [200, 201]])
        self.assertEqual(model["lengths"]["cds"], 6)

    def test_stop_codon_already_inside_cds_is_not_double_counted(self):
        tx = record()
        tx["cds"] = [[105, 110], [200, 209]]
        self.assertEqual(segment_transcript(tx)["lengths"]["cds"], 14)

    def test_cds_with_unannotated_internal_gap_is_rejected(self):
        tx = record()
        tx["cds"] = [[105, 109], [200, 206]]
        with self.assertRaisesRegex(ValueError, "discontinuous_spliced_cds"):
            segment_transcript(tx)

    def test_intronic_coding_features_rejected(self):
        tx = record()
        tx["cds"] = [[105, 200]]
        with self.assertRaisesRegex(ValueError, "coding_feature_outside_exons"):
            segment_transcript(tx)

    def test_missing_utr_and_partial_cds_are_reported(self):
        tx = record()
        tx["cds"] = [[100, 110], [200, 206]]
        with self.assertRaisesRegex(ValueError, "missing_or_short_utr_or_cds"):
            segment_transcript(tx)
        tx["partial_cds"] = True
        with self.assertRaisesRegex(ValueError, "partial_cds_annotation"):
            segment_transcript(tx)

    def test_longest_mature_transcript_not_longest_genomic_span_or_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.gtf"
            lines = []
            # Long genomic span/long CDS has only 20 nt of mature RNA. The
            # correct representative has 30 nt of RNA and just 6 nt of CDS.
            for tx, features in {
                "long_genomic": [("exon", 101, 110), ("exon", 1001, 1010),
                                 ("CDS", 105, 110), ("CDS", 1001, 1008)],
                "long_mature": [("exon", 201, 230), ("CDS", 211, 216)],
                "no_utr": [("exon", 301, 400), ("CDS", 301, 400)],
            }.items():
                for feature, start, end in features:
                    lines.append(f'chrTest\ttest\t{feature}\t{start}\t{end}\t.\t+\t.\tgene_id "g"; transcript_id "{tx}";\n')
            # Include a second excluded gene to verify exclusion accounting.
            lines.append('chrTest\ttest\texon\t501\t510\t.\t+\t.\tgene_id "noncoding"; transcript_id "nc";\n')
            path.write_text("".join(lines[::-1]))
            models, rows, summary = prepare_annotation(path)
            self.assertEqual(models[0]["transcript_id"], "long_mature")
            self.assertEqual(summary["genes_included"], 1)
            self.assertEqual(summary["genes_excluded"], 1)
            self.assertEqual(summary["transcript_exclusions"]["noncoding"], 1)
            self.assertEqual(summary["transcript_exclusions"]["missing_or_short_utr_or_cds"], 1)


@unittest.skipIf(np is None, "numpy required for signal integration tests")
class SignalTests(unittest.TestCase):
    def test_fractional_bin_mass_is_conserved(self):
        bw = FakeBigWig([(10, 11, 2.0), (12, 13, 4.0)])
        profile = scaled_region_profile(bw, "chrTest", [[10, 13]], "+", 2)
        np.testing.assert_allclose(profile, [2 / 1.5, 4 / 1.5])
        self.assertAlmostEqual(sum(profile) * 1.5, 6)

    def test_short_region_more_bins_than_bases(self):
        bw = FakeBigWig([(10, 11, 7.0)])
        np.testing.assert_allclose(scaled_region_profile(bw, "chrTest", [[10, 11]], "+", 100), 7)

    def test_minus_strand_reverses_and_skips_intron(self):
        bw = FakeBigWig([(10, 11, 1.0), (11, 12, 2.0), (12, 20, 999.0),
                         (20, 21, 3.0), (21, 22, 4.0)])
        profile = scaled_region_profile(bw, "chrTest", [[20, 22], [10, 12]], "-", 4)
        np.testing.assert_allclose(profile, [4, 3, 2, 1])

    def test_empty_and_partial_bigwig_coverage(self):
        bw = FakeBigWig([(0, 12, 5.0)])
        np.testing.assert_allclose(scaled_region_profile(bw, "chrTest", [[10, 14]], "+", 4), [5, 5, 0, 0])

    def test_actual_bigwig_to_density_and_gene_signal(self):
        try:
            import pyBigWig
            import matplotlib
        except ImportError:
            self.skipTest("pyBigWig/matplotlib required for end-to-end metaplot test")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tx = record("-")
            tx.update(exons=[[100, 112], [200, 210]], cds=[[106, 112], [200, 205]], stop_codons=[[103, 106]])
            model = segment_transcript(tx)
            annotation = root / "models.jsonl"
            annotation.write_text(json.dumps(model) + "\n")
            summary = root / "summary.json"
            summary.write_text(json.dumps({"region_median_lengths": model["lengths"], "selection_policy": "test"}))
            # Gene is negative strand, so positive track's huge decoy signal
            # must never appear. The minus RNA 5' endpoint is in the first UTR5 bin.
            paths = []
            for strand, start, value in [("plus", 207, 999.0), ("minus", 209, 2.0)]:
                path = root / f"{strand}.bw"
                with pyBigWig.open(str(path), "w") as bw:
                    bw.addHeader([("chrTest", 300)])
                    bw.addEntries(["chrTest"], [start], ends=[start + 1], values=[value])
                paths.append(str(path))
            output = Named(plot=str(root / "plot.pdf"), data=str(root / "profile.tsv"),
                           stats=str(root / "stats.json"), gene_signal=str(root / "genes.tsv"))
            smk = SimpleNamespace(scriptdir=str(SCRIPTS), wildcards=SimpleNamespace(sample="synthetic"),
                input=Named(bws=paths, models=str(annotation), annotation_summary=str(summary)), output=output,
                params=SimpleNamespace(bins={"utr5": 5, "cds": 14, "utr3": 3},
                    tracks=["crosslinkCPM.noblacklist"], track_specs=[{"plus": 0, "minus": 1}],
                    axis_scale="median_length", axis_widths=[1, 1, 1]))
            runpy.run_path(str(SCRIPTS / "compute_segmented_metaplot.py"), init_globals={"snakemake": smk})
            with open(output.data) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertAlmostEqual(sum(float(r["probability_density"]) * float(r["x_width"]) for r in rows), 1)
            self.assertAlmostEqual(float(rows[0]["mean_cpm"]), 2)
            self.assertTrue(all(float(r["mean_cpm"]) == 0 for r in rows[1:]))
            with open(output.gene_signal) as handle:
                gene = next(csv.DictReader(handle, delimiter="\t"))
            self.assertAlmostEqual(float(gene["signal_mass"]), 2)
            self.assertGreater(Path(output.plot).stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
