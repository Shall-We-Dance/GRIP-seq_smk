"""Regression tests for validated, independent BAM-coverage reuse."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pysam
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('bam_coverage', ROOT / 'workflow/rules/scripts/bam_coverage.py')
coverage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coverage)


class NamedInputs(dict):
    __getattr__ = dict.__getitem__


class BamCoverageTests(unittest.TestCase):
    def fixture(self, directory, chromosome='chrTest', length=100):
        directory = Path(directory)
        bam = directory / 'aligned reads.bam'
        with pysam.AlignmentFile(str(bam), 'wb', header={
            'HD': {'VN': '1.6', 'SO': 'coordinate'},
            'SQ': [{'SN': 'chrTest', 'LN': 100}],
        }):
            pass
        source = directory / 'baseline coverage.bw'
        with pyBigWig.open(str(source), 'w') as track:
            track.addHeader([(chromosome, length)])
            track.addEntries([chromosome], [10], ends=[12], values=[5.0])
        return SimpleNamespace(
            input=NamedInputs(bam=str(bam), reuse=[str(source)]),
            output=NamedInputs(bw=str(directory / 'new results' / 'coverage.bw')),
            log=[str(directory / 'new results' / 'coverage.log')],
            threads=1,
        )

    def test_copy_preserves_bytes_but_not_source_inode(self):
        with tempfile.TemporaryDirectory() as temporary:
            sm = self.fixture(temporary)
            source = Path(sm.input.reuse[0])
            original = source.read_bytes()
            with patch.object(coverage.subprocess, 'run') as run:
                coverage.create_coverage(sm)
                run.assert_not_called()
            output = Path(sm.output.bw)
            self.assertEqual(output.read_bytes(), original)
            self.assertNotEqual(output.stat().st_ino, source.stat().st_ino)
            self.assertFalse(output.is_symlink())
            self.assertEqual(json.loads(Path(sm.log[0]).read_text())['mode'], 'reused_validated_bam_coverage')
            output.write_bytes(b'replaced output')
            self.assertEqual(source.read_bytes(), original)

    def test_mismatched_chromosome_names_or_lengths_are_rejected(self):
        for chromosome, length in [('chrDifferent', 100), ('chrTest', 101)]:
            with self.subTest(chromosome=chromosome, length=length), tempfile.TemporaryDirectory() as temporary:
                sm = self.fixture(temporary, chromosome, length)
                with self.assertRaisesRegex(ValueError, 'chromosome names and lengths'):
                    coverage.create_coverage(sm)
                self.assertFalse(Path(sm.output.bw).exists())

    def test_fallback_keeps_original_coverage_parameters(self):
        for blacklist in (False, True):
            with self.subTest(blacklist=blacklist), tempfile.TemporaryDirectory() as temporary:
                sm = self.fixture(temporary)
                sm.input['reuse'] = []
                sm.threads = 3
                if blacklist:
                    sm.input['bl'] = str(Path(temporary) / 'blacklist intervals.bed')
                expected = ['bamCoverage', '-b', sm.input.bam, '-o', sm.output.bw,
                            '--binSize', '1', '--normalizeUsing', 'CPM',
                            '--numberOfProcessors', '3']
                if blacklist:
                    expected += ['--blackListFileName', sm.input.bl]
                with patch.object(coverage.subprocess, 'run') as run:
                    coverage.create_coverage(sm)
                self.assertEqual(run.call_args.args[0], expected)
                self.assertTrue(run.call_args.kwargs['check'])
                self.assertEqual(run.call_args.kwargs['stderr'], coverage.subprocess.STDOUT)


if __name__ == '__main__':
    unittest.main()
