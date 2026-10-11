"""Regression: retain significant de novo results in the complete sample table."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow/rules/scripts"
sys.path.insert(0, str(SCRIPTS))
from summarize_motif_offsets import read_results, sample_summary_rows


class MotifSummaryTests(unittest.TestCase):
    def fixture(self, root, sample, consensus, qvalue=.01):
        motifs = [
            {"name": "DRACH", "sequence": "DRACH", "anchor_pos": 3, "source": "predefined",
             "effective_site_pairs": 100, "supporting_reads": 1000, "status": "completed",
             "peak_frequency_offsets_nt": [0], "peak_excess_over_background_offsets_nt": [0],
             "significant_enrichment_offsets_q_lt_0p05": []},
            {"name": "MEME_shared_name", "sequence": consensus, "anchor_pos": 2,
             "source": "de_novo_consensus_holdout", "effective_site_pairs": 20,
             "supporting_reads": 200, "status": "completed", "peak_frequency_offsets_nt": [-1],
             "peak_excess_over_background_offsets_nt": [-1],
             "significant_enrichment_offsets_q_lt_0p05": [-1]},
        ]
        statistics = root / f"{sample}.json"
        statistics.write_text(json.dumps({"sample": sample, "motifs": motifs,
                                         "multiple_testing": "sample-level family"}))
        table = root / f"{sample}.tsv"
        fields = ["motif", "iupac", "anchor_pos_1based", "source", "effective_site_pairs",
                  "coordinate_population", "motif_anchor_offset_from_crosslink_rna",
                  "foreground_site_fraction", "background_site_fraction", "qvalue_bh_sample"]
        with table.open("w") as out:
            writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for motif in motifs:
                for offset in [-1, 0, 1]:
                    denovo_peak = motif["source"] != "predefined" and offset == -1
                    writer.writerow(dict(zip(fields, [motif["name"], motif["sequence"], motif["anchor_pos"],
                        motif["source"], motif["effective_site_pairs"], "all_eligible_crosslinks", offset,
                        .5 if denovo_peak else .1, .05, qvalue if denovo_peak else 1.])))
        return table, statistics

    def test_summary_keeps_sample_specific_de_novo_while_plot_definitions_are_predefined(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = self.fixture(root, "sample_A", "TGAC", .003)
            b = self.fixture(root, "sample_B", "CGAT", .02)
            samples = ["sample_A", "sample_B"]
            profiles, reports, definitions, offsets = read_results(samples, [a[0], b[0]], [a[1], b[1]])
            self.assertEqual(definitions, {"DRACH": ("DRACH", 3)})
            self.assertEqual(offsets, {"DRACH": [-1, 0, 1]})
            rows = sample_summary_rows(samples, reports, profiles, {"comparison": {"IP": samples}})
            self.assertEqual(len(rows), 4)
            denovo = [row for row in rows if row["source"] == "de_novo_consensus_holdout"]
            self.assertEqual([row["iupac"] for row in denovo], ["TGAC", "CGAT"])
            self.assertEqual([row["minimum_qvalue_bh_sample"] for row in denovo], [.003, .02])
            self.assertEqual([json.loads(row["significant_enrichment_offsets_q_lt_0p05"]) for row in denovo], [[-1], [-1]])
            self.assertEqual([row["effective_site_pairs"] for row in denovo], [20, 20])


if __name__ == "__main__":
    unittest.main()
