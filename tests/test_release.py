"""Public exports must preserve public sources while excluding private run data."""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('build_release', ROOT / 'tools/build_release.py')
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class ReleaseTests(unittest.TestCase):
    def source(self, directory):
        root = Path(directory) / 'source'
        root.mkdir()
        for name in release.TOP:
            (root / name).write_text('Public source\n')
        (root / 'LICENSE').write_text('Copyright (c) 2026 Public Author\n')
        (root / 'README.md').write_text('Source: https://github.com/public/gripseq\n')
        (root / '.gitignore').write_text('/data/\n/config.yaml\n')
        for relative in ['workflow/Snakefile', 'workflow/rules/envs/test.yaml', 'tests/test_public.py', 'tools/public.py', '.github/workflows/ci.yaml']:
            p = root / relative; p.parent.mkdir(parents=True, exist_ok=True); p.write_text('# public source\n')
        return root

    def export(self, directory, root):
        with contextlib.redirect_stdout(io.StringIO()):
            return release.export(directory, root=root)

    def test_allowlist_archive_manifest_and_license(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self.source(temporary)
            for relative in ['.git/config', 'config.yaml', 'results/metadata.json', 'workflow/config.yaml', 'workflow/results/private.py', 'workflow/rules/__pycache__/private.py', 'workflow/raw.bam', 'tests/provenance.json']:
                p = root / relative; p.parent.mkdir(parents=True, exist_ok=True); p.write_text('samples: {}\n' if p.name == 'config.yaml' else 'private metadata\n')
            destination = Path(temporary) / 'public.v2'
            result = self.export(destination, root)
            self.assertEqual(result['archive'], str(destination) + '.tar.gz')
            self.assertEqual((destination / 'LICENSE').read_bytes(), (root / 'LICENSE').read_bytes())
            manifest = json.loads((destination / 'RELEASE_MANIFEST.json').read_text())['files']
            actual = {str(p.relative_to(destination)) for p in destination.rglob('*') if p.is_file()}
            self.assertEqual(actual, set(manifest) | {'RELEASE_MANIFEST.json'})
            for name, digest in manifest.items():
                self.assertEqual(hashlib.sha256((destination / name).read_bytes()).hexdigest(), digest)
                self.assertFalse(release.excluded(Path(name)))
            with tarfile.open(result['archive']) as tar:
                archived = {str(Path(m.name).relative_to(destination.name)) for m in tar.getmembers() if m.isfile()}
                self.assertEqual(archived, actual)
                self.assertTrue(all(not m.issym() and not m.islnk() for m in tar.getmembers()))
                for member in tar.getmembers():
                    self.assertEqual((member.uid, member.gid, member.uname, member.gname, member.mtime), (0, 0, '', '', 0))

    def test_private_content_fails_before_creating_output(self):
        # Assemble fixtures so the public test source contains no private path
        # or credential literal that could match the exporter's own checks.
        fixtures = ['/' + 'home' + '/private_user/raw.fq', 'CKDL' + '123456789', 'api_' + 'key = ' + chr(34) + 'example-secret' + chr(34)]
        for text in fixtures:
            with self.subTest(text_kind=text.split('/')[0]), tempfile.TemporaryDirectory() as temporary:
                root = self.source(temporary)
                (root / 'workflow/Snakefile').write_text(text)
                destination = Path(temporary) / 'public'
                with self.assertRaisesRegex(ValueError, 'Potential private content'):
                    self.export(destination, root)
                self.assertFalse(destination.exists())
                self.assertFalse(Path(str(destination) + '.tar.gz').exists())

    def test_local_sample_names_are_checked_without_exporting_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self.source(temporary)
            sample = 'local_' + 'donor_007'
            (root / 'config.yaml').write_text('samples:\n  ' + sample + ': {}\n')
            (root / 'workflow/Snakefile').write_text('# ' + sample + '\n')
            destination = Path(temporary) / 'public'
            with self.assertRaisesRegex(ValueError, 'Potential private content'):
                self.export(destination, root)
            self.assertFalse(destination.exists())

    def test_existing_archive_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self.source(temporary)
            destination = Path(temporary) / 'public'
            archive = Path(str(destination) + '.tar.gz')
            archive.write_bytes(b'previous release')
            with self.assertRaisesRegex(ValueError, 'never overwritten'):
                self.export(destination, root)
            self.assertEqual(archive.read_bytes(), b'previous release')
            self.assertFalse(destination.exists())

    def test_symlink_sources_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self.source(temporary)
            target = Path(temporary) / 'private.txt'; target.write_text('outside source')
            (root / 'workflow/linked.py').symlink_to(target)
            destination = Path(temporary) / 'public'
            with self.assertRaisesRegex(ValueError, 'Symlink'):
                self.export(destination, root)
            self.assertFalse(destination.exists())

    def test_exporter_does_not_flag_its_own_code(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self.source(temporary)
            (root / 'tools/build_release.py').write_bytes(Path(release.__file__).read_bytes())
            (root / 'tests/test_release.py').write_bytes(Path(__file__).read_bytes())
            result = self.export(Path(temporary) / 'public', root)
            self.assertTrue(Path(result['archive']).is_file())


if __name__ == '__main__':
    unittest.main()
