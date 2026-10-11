"""Analytic candidate coordinates, Input normalization, set logic and site density."""
import gzip
import json
from pathlib import Path
import sys
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

import pysam
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"workflow/rules/scripts"))
import m6a_sites as m6a
import compare_m6a_sites as compare


class Named(dict):
    __getattr__ = dict.__getitem__


class CandidateTests(unittest.TestCase):
    def reference(self, directory):
        sequence = list("C"*300)
        sequence[41:46] = "GGACT"     # RNA A at 43; plus crosslink at 40
        sequence[135:140] = "AGTCC"  # reverse complement GGACT: A at 137
        path = Path(directory)/"reference.fa"
        path.write_text(">chrTest\n"+"".join(sequence)+"\n")
        pysam.faidx(str(path))
        return path

    def data(self):
        keys = [("chrTest", 40, "+"), ("chrTest", 140, "-")]
        sites = {key: dict(count=12, left=2, right=3, r2_positions={key[1]+(1 if key[2] == "+" else -1)}) for key in keys}
        return keys, sites

    def test_both_strands_input_normalization_and_unselected_offsets(self):
        with tempfile.TemporaryDirectory() as tmp, pysam.FastaFile(str(self.reference(tmp))) as fasta:
            keys, sites = self.data()
            rows, offsets, eligible = m6a.infer_sites(sites, dict.fromkeys(keys, 12), {}, fasta, "example", 100, 100, {})
            self.assertEqual([(row["pos0"], row["RNA_5mer"], row["passes_m6A"]) for row in rows], [(43, "GGACT", 1), (137, "GGACT", 1)])
            self.assertEqual(offsets[-3], 2)
            self.assertEqual(eligible, 2)
            self.assertAlmostEqual(rows[0]["IP_Input_fold"], 25)
            # Input filtering changes calls, never the independent offset curve.
            failed, offsets2, _ = m6a.infer_sites(sites, dict.fromkeys(keys, 12), dict.fromkeys(keys, 12), fasta, "example", 100, 100, {})
            self.assertFalse(any(row["passes_m6A"] for row in failed))
            self.assertEqual(offsets, offsets2)

    def test_read_threshold_is_strict_and_zero_depth_is_not_estimable(self):
        with tempfile.TemporaryDirectory() as tmp, pysam.FastaFile(str(self.reference(tmp))) as fasta:
            keys, sites = self.data()
            for value, expected in [(10, 0), (11, 1)]:
                for site in sites.values():
                    site["count"] = value
                rows, _, _ = m6a.infer_sites(sites, dict.fromkeys(keys, value), {}, fasta, "example", 100, 100, {})
                self.assertEqual([row["passes_m6A"] for row in rows], [expected]*2)
            rows, _, _ = m6a.infer_sites(sites, dict.fromkeys(keys, 11), {}, fasta, "example", 100, 0, {})
            self.assertTrue(all(row["IP_Input_fold"] is None and not row["passes_m6A"] for row in rows))

    def test_pseudocount_is_applied_before_library_depth_scaling(self):
        self.assertAlmostEqual(m6a.enrichment(12, 0, 200, 100), 12.5)
        self.assertIsNone(m6a.enrichment(0, 0, 0, 100))

    def test_atlas_bed_validates_strand_base_and_half_open_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp, pysam.FastaFile(str(self.reference(tmp))) as fasta:
            path = Path(tmp)/"atlas.bed"
            path.write_text("chrTest\t43\t44\ta\t0\t+\nchrTest\t137\t138\tb\t0\t-\n")
            self.assertEqual(m6a.load_atlas(path, fasta), {("chrTest", 43, "+"), ("chrTest", 137, "-")})
            for invalid in ["chrTest\t44\t45\ta\t0\t+\n", "chrTest\t43\t43\ta\t0\t+\n", "chrTest\t43\t44\ta\t0\t-\n"]:
                path.write_text(invalid)
                with self.assertRaises(ValueError):
                    m6a.load_atlas(path, fasta)

    def test_tolerance_matching_is_one_to_one_and_strand_specific(self):
        first = {("chrTest", 10, "+"), ("chrTest", 12, "+"), ("chrTest", 30, "-")}
        second = {("chrTest", 11, "+"), ("chrTest", 30, "+")}
        self.assertEqual(m6a.match_count(first, second, 1), 1)
        self.assertEqual(m6a.match_count(first, second, 0), 0)
        second.add(("chrTest", 13, "+"))
        self.assertEqual(m6a.match_count(first, second, 1), 2)

    def test_replicate_consensus_and_disjoint_memberships(self):
        a, b, c = [("chrTest", p, "+") for p in (10, 20, 30)]
        self.assertEqual(m6a.consensus([{a,b}, {b,c}], 2), {b})
        labels, patterns = compare.memberships({"A": {a,b}, "B": {b,c}})
        self.assertEqual(sum(patterns.values()), 3)
        self.assertEqual(patterns[(1,1)], 1)

    def test_fractional_gene_assignment_conserves_site_weight_and_density(self):
        model = dict(chrom="chrTest", strand="+", exons=[[0,30],[100,130]],
                     lengths=dict(utr5=10, cds=40, utr3=10), gene_id="g1", transcript_id="t1")
        model2 = dict(model, gene_id="g2", transcript_id="t2")
        ann = m6a.SiteAnnotation([model,model2])
        sites = {("chrTest", 5, "+"), ("chrTest", 125, "+"), ("chrTest", 70, "+")}
        counts, density, assigned = ann.profile(sites)
        self.assertEqual(assigned, 2)
        self.assertAlmostEqual(sum(counts), 2)
        self.assertAlmostEqual(sum(density)/30, 1)
        self.assertEqual(ann.profile(set())[1], [None]*90)
        model.update(strand="-")
        self.assertAlmostEqual(m6a.SiteAnnotation([model]).annotate(("chrTest", 125, "-"))[0][2], 0.4)

    def test_comparison_exports_empty_sets_atlas_and_insufficient_replication(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = self.reference(tmp)
            calls, offsets, stats = [], [], []
            for name, sites in [("ip1", [43]), ("ip2", [])]:
                call = tmp/f"{name}.tsv"
                m6a.write_tsv(call, [dict(chrom="chrTest", pos0=p, rna_strand="+") for p in sites], ["chrom", "pos0", "rna_strand"])
                off = tmp/f"{name}.offset.tsv"
                m6a.write_tsv(off, [dict(crosslink_minus_A_RNA_nt=-3, site_anchor_pairs=1, eligible_RT_sites=1, fraction=1)],
                             ["crosslink_minus_A_RNA_nt", "site_anchor_pairs", "eligible_RT_sites", "fraction"])
                stat = tmp/f"{name}.json";stat.write_text("{}")
                calls.append(str(call));offsets.append(str(off));stats.append(str(stat))
            atlas = tmp/"atlas.bed";atlas.write_text("chrTest\t43\t44\ta\t0\t+\n")
            gtf = tmp/"test.gtf"
            gtf.write_text('chrTest\ttest\texon\t1\t300\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
                           'chrTest\ttest\tCDS\t51\t250\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
            sm = SimpleNamespace(input=Named(calls=calls, offsets=offsets, stats=stats, fasta=str(fasta), atlas=[str(atlas)], gtf=[str(gtf)]),
                                 output=Named(report=str(tmp/"report")),
                                 params=Named(samples=["ip1", "ip2"], settings=dict(groups={"condition": ["ip1", "ip2"], "single": ["ip1"]})))
            compare.run(sm)
            report = json.loads((tmp/"report/report.json").read_text())
            single = next(row for row in report["group_definitions"] if row["scheme"] == "replicated" and row["group"] == "single")
            self.assertEqual(single["status"], "not_estimable_insufficient_libraries")
            self.assertIsNone(single["n_sites"])
            profiles = m6a.read_tsv(tmp/"report/site_metagene.tsv")
            density = [float(row["density"]) for row in profiles if row["scheme"] == "union" and row["label"] == "condition" and row["subset"] == "all"]
            self.assertAlmostEqual(sum(density)/30, 1)
            self.assertTrue((tmp/"report/UpSet.Atlas.union.pdf").is_file())
            self.assertEqual(m6a.read_tsv(tmp/"report/Atlas_unmatched_within5.tsv"), [])

    def test_snakemake_executes_calling_and_comparison_with_synthetic_signals(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            shutil.copytree(ROOT/"workflow", tmp/"workflow", ignore=shutil.ignore_patterns("__pycache__"))
            fasta = self.reference(tmp)
            config = yaml.safe_load((ROOT/"config.example.yaml").read_text())
            config["reference"]["fasta"] = str(fasta)
            config["reference"]["gtf"] = None
            config["filter_blacklist"] = False
            config["samples"] = {"ip": {"bam": "unused-ip.bam"}, "input": {"bam": "unused-input.bam"}}
            for module in ("metaplot", "de_novo_motif", "motif_distribution", "bigwigCompare"):
                config[module]["enable"] = False
            config["bigwigCompare"]["groups"] = {}
            config["m6a"].update(enable=True, controls={"ip": "input"}, groups={"example": ["ip"]})
            (tmp/"config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
            for sample in ("ip", "input"):
                directory = tmp/f"results/signal/{sample}"
                directory.mkdir(parents=True)
                (directory/f"{sample}.noblacklist.signal.stats.json").write_text(json.dumps(dict(accepted_read2=100)))
                with gzip.open(directory/f"{sample}.crosslink.noblacklist.counts.bed.gz", "wt") as handle:
                    handle.write("chrTest\t40\t41\tx\t12\t+\nchrTest\t140\t141\ty\t12\t-\n" if sample == "ip" else "chrTest\t200\t201\tother\t100\t+\n")
            directory = tmp/"results/signal/ip"
            with gzip.open(directory/"ip.R2firstbase.noblacklist.counts.bed.gz", "wt") as handle:
                handle.write("chrTest\t41\t42\tx\t12\t+\nchrTest\t139\t140\ty\t12\t-\n")
            with gzip.open(directory/"ip.R2_to_crosslink.noblacklist.counts.tsv.gz", "wt") as handle:
                handle.write("chrom\tr2_pos0\tcrosslink_pos0\trna_strand\tcount\nchrTest\t41\t40\t+\t12\nchrTest\t139\t140\t-\t12\n")
            result = subprocess.run([shutil.which("snakemake"), "-s", "workflow/Snakefile", "results/m6a/comparison",
                                     "--cores", "1", "--nocolor", "--allowed-rules", "call_m6a_candidates", "compare_m6a_candidates"],
                                    cwd=tmp, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
            self.assertEqual(result.returncode, 0, result.stdout)
            rows = m6a.read_tsv(tmp/"results/m6a/ip/candidates.tsv")
            self.assertEqual({m6a.site_key(row) for row in rows}, {("chrTest",43,"+"), ("chrTest",137,"-")})
            self.assertTrue((tmp/"results/m6a/comparison/UpSet.union.pdf").is_file())


if __name__ == "__main__":
    unittest.main()
