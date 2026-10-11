"""Compare the vectorized writer to real deepTools on synthetic gapped reads."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pysam
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]


class CoverageEquivalenceTests(unittest.TestCase):
    def test_basewise_cpm_matches_deeptools_with_and_without_blacklist(self):
        executable = shutil.which("bamCoverage")
        self.assertIsNotNone(executable, "Install deepTools to run coverage equivalence tests")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bam = tmp/"reads.bam"
            with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                      "SQ": [{"SN": "chrTest", "LN": 1000}, {"SN": "chrEmpty", "LN": 100}]}) as handle:
                for i, (position, cigar) in enumerate([(10, "20M"), (10, "20M"), (15, "10M5N10M"),
                                                       (80, "5S20M"), (90, "10M2D10M"), (200, "20M")]):
                    read = pysam.AlignedSegment()
                    read.query_name=f"synthetic{i}"; read.reference_id=0; read.reference_start=position
                    read.flag=0; read.mapping_quality=60; read.cigarstring=cigar
                    read.query_sequence="A"*read.infer_query_length();handle.write(read)
            pysam.index(str(bam))
            blacklist = tmp/"blacklist.bed";blacklist.write_text("chrTest\t195\t240\n")
            python = str(Path(executable).resolve().parent/"python")
            if not Path(python).exists():
                python = sys.executable
            for blocked in (False, True):
                with self.subTest(blacklist=blocked):
                    first, second = tmp/"original.bw", tmp/"vectorized.bw"
                    common = ["-b", str(bam), "--binSize", "1", "--normalizeUsing", "CPM", "--numberOfProcessors", "1"]
                    if blocked:
                        common += ["--blackListFileName", str(blacklist)]
                    for runner, output in [([executable], first), ([python, str(ROOT/"workflow/rules/scripts/fast_bam_coverage.py")], second)]:
                        result = subprocess.run(runner+common+["-o", str(output)], capture_output=True, text=True, timeout=60)
                        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                    with pyBigWig.open(str(first)) as a, pyBigWig.open(str(second)) as b:
                        self.assertEqual(a.chroms(), b.chroms())
                        for chrom, length in a.chroms().items():
                            np.testing.assert_allclose(a.values(chrom,0,length), b.values(chrom,0,length), equal_nan=True, rtol=0, atol=0)


if __name__ == "__main__":
    unittest.main()
