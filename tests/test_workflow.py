"""Offline DAG smoke tests; no sequencing tools or real data are executed."""

import ast
import copy
import gzip
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import yaml


REPO = Path(__file__).resolve().parents[1]


class WorkflowSmokeTest(unittest.TestCase):
    def test_python_and_yaml_syntax(self):
        for path in (REPO / "workflow").rglob("*.py"):
            ast.parse(path.read_text(), filename=str(path))
        for path in [REPO / "config.example.yaml", REPO / "environment.yaml",
                     *(REPO / "workflow/rules/envs").glob("*.yaml")]:
            self.assertIsInstance(yaml.safe_load(path.read_text()), dict, str(path))

    def test_dag_with_optional_modules(self):
        executable = shutil.which("snakemake")
        self.assertIsNotNone(executable, "Install environment.yaml to run DAG tests")
        template = yaml.safe_load((REPO / "config.example.yaml").read_text())
        for variant in ("all", "core", "no_compare", "no_blacklist", "bam_input", "m6a_only", "no_cleanup", "keep_raw"):
            with self.subTest(variant=variant), tempfile.TemporaryDirectory(prefix="gripseq-test-") as tmp:
                root = Path(tmp)
                # Use only distributable code: the real config/results are never read.
                shutil.copytree(REPO / "workflow", root / "workflow")
                config = copy.deepcopy(template)
                if variant == "all":
                    config["motif_anchoring"]["enable"] = True
                    config["m6a"]["enable"] = True
                if variant == "m6a_only":
                    config["m6a"]["enable"] = True
                    for module in ("de_novo_motif", "motif_distribution", "metaplot", "bigwigCompare"):
                        config[module]["enable"] = False
                if variant == "no_cleanup":
                    config["adapter_cleanup"]["enable"] = False
                if variant == "keep_raw":
                    config["preprocessing"]["keep_merged_fastq"] = True
                if variant == "core":
                    for module in ("bigwigCompare", "motif_anchoring", "de_novo_motif", "motif_distribution", "metaplot"):
                        config[module]["enable"] = False
                    config["filter_blacklist"] = False
                if variant == "no_compare":
                    # Keep groups populated to catch disabled-module target regressions.
                    config["bigwigCompare"]["enable"] = False
                if variant == "no_blacklist":
                    config["filter_blacklist"] = False
                config["blacklist"].pop("url", None)
                config["blacklist"]["path"] = "resources/test.blacklist.bed"
                (root / "resources").mkdir()
                (root / "resources/test.blacklist.bed").write_text("chrTest\t0\t1\n")
                for sample in config["samples"].values():
                    for paths in sample.values():
                        for filename in paths:
                            path = root / filename
                            path.parent.mkdir(parents=True, exist_ok=True)
                            with gzip.open(path, "wt") as handle:
                                handle.write("@synthetic\nACGT\n+\nIIII\n")
                if variant == "bam_input":
                    for name in config["samples"]:
                        path = root / "data" / f"{name}.bam"
                        path.touch()
                        config["samples"][name] = {"bam": str(path)}
                reference = root / config["reference"]["fasta"]
                reference.parent.mkdir(parents=True, exist_ok=True)
                reference.write_text(">chrTest\n" + "ACGT" * 25 + "\n")
                (root / config["reference"]["gtf"]).write_text(
                    'chrTest\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "test";\n'
                )
                (root / config["reference"]["star_index"]).mkdir()
                (root / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
                result = subprocess.run(
                    [executable, "-s", "workflow/Snakefile", "--dry-run", "--cores", "2",
                     "--printshellcmds", "--nocolor"],
                    cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    timeout=120,
                )
                self.assertEqual(result.returncode, 0, result.stdout)
                self.assertIn("rule all:", result.stdout)
                for rule in ("rule star_align_unique:", "rule multiqc:"):
                    if variant == "bam_input":
                        self.assertNotIn(rule, result.stdout)
                    else:
                        self.assertIn(rule, result.stdout)
                self.assertNotIn("cleanup_tmp", result.stdout)
                for enabled, rule in (
                    (config["bigwigCompare"]["enable"], "rule bigwig_compare_log2:"),
                    (config["filter_blacklist"], "rule bam_to_bigwig_cpm_blacklist:"),
                    (config["motif_anchoring"]["enable"], "rule motif_anchor:"),
                    (config["de_novo_motif"]["enable"], "rule meme_de_novo_pool:"),
                    (config["metaplot"]["enable"], "rule plot_metaplot_integrated:"),
                    (config["motif_distribution"]["enable"], "rule motif_offset_distribution:"),
                    (config["m6a"]["enable"], "rule call_m6a_candidates:"),
                    (config["m6a"]["enable"], "rule compare_m6a_candidates:"),
                    (config["adapter_cleanup"]["enable"] and variant != "bam_input", "rule adapter_cleanup:"),
                ):
                    if enabled:
                        self.assertIn(rule, result.stdout)
                    else:
                        self.assertNotIn(rule, result.stdout)


if __name__ == "__main__":
    unittest.main()
