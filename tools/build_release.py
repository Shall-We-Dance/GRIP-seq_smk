"""Export a reviewed source allowlist without local configuration or run data."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
TOP = ['README.md', 'METHODS.md', 'LICENSE', 'CHANGELOG.md', 'config.example.yaml', 'environment.yaml', '.gitignore']
TREES = ['workflow', 'tests', 'tools', '.github']
SUFFIXES = {'.py', '.smk', '.yaml', '.yml', '.md', '.sh'}
EXCLUDED_PARTS = {'.git', '.snakemake', '.codex', '.agents', '__pycache__', '.pytest_cache',
                  'logs', 'results', 'diagnostics', 'data', '01.RawData', 'dist'}
# Case matters for private storage roots: a portable /data/ ignore rule is safe.
# Require a path boundary so a public URL containing the same path segment passes.
PRIVATE_PATH = re.compile(r'(?<![A-Za-z0-9:/.-])/(?:home|Users|FAST|DATA|BIG|public)/[^\s"\x27]+')
DELIVERY_ID = re.compile(r'\bCKDL[0-9]{6,}[A-Za-z0-9_-]*\b', re.I)
CREDENTIAL = re.compile(r'\b(?:api[_-]?key|password|access[_-]?token|secret[_-]?key)["\x27]?\s*[:=]\s*["\x27][^"\x27\n]+["\x27]', re.I)


def excluded(relative):
    if any(part in EXCLUDED_PARTS for part in relative.parts):
        return True
    name = relative.name
    return (name == 'config.yaml' or (name.startswith('config.') and name.endswith('.local.yaml'))
            or name.startswith('.env') or name in {'metadata.json', 'provenance.json'})


def source_files(root):
    """Select source-like files, rejecting links instead of following their targets."""
    files = []
    for name in TOP:
        candidate = root / name
        if not candidate.is_file():
            raise ValueError(f'Required public source is missing: {name}')
        files.append(candidate)
    for tree in TREES:
        base = root / tree
        if base.is_symlink():
            raise ValueError(f'Symlink source directories are not exported: {tree}')
        for candidate in base.rglob('*'):
            relative = candidate.relative_to(root)
            if excluded(relative):
                continue
            if candidate.is_symlink():
                raise ValueError(f'Symlink source files are not exported: {relative}')
            if candidate.is_file() and (candidate.suffix in SUFFIXES or candidate.name == 'Snakefile'):
                files.append(candidate)
    return sorted(set(files))


def private_sample_names(root):
    """Read labels from excluded local configs; never embed an experiment's names.

    PyYAML is already an explicit workflow-environment dependency. No local
    config means this utility can still run in a clean public checkout without it.
    """
    configs = ([root / 'config.yaml'] if (root / 'config.yaml').exists() else [])
    configs += sorted(root.glob('config.*.local.yaml'))
    if not configs:
        return set()
    try:
        import yaml
    except ImportError as exc:
        raise ValueError('Run in the workflow environment (PyYAML is needed to inspect local sample names)') from exc
    public_config = yaml.safe_load((root / 'config.example.yaml').read_text()) or {}
    public_names = set(public_config.get('samples', {})) if isinstance(public_config, dict) else set()
    names = set()
    for path in configs:
        local = yaml.safe_load(path.read_text()) or {}
        if not isinstance(local, dict):
            raise ValueError(f'Expected a mapping in local config: {path.name}')
        names.update(str(name) for name in local.get('samples', {}))
    return names - public_names


def inspect_sources(root):
    """Read and validate every byte before creating a release directory."""
    snapshots = {}
    private_names = private_sample_names(root)
    for source in source_files(root):
        relative = source.relative_to(root)
        if source.is_symlink():
            raise ValueError(f'Symlink source files are not exported: {relative}')
        content = source.read_bytes()
        try:
            text = content.decode('utf-8')
        except UnicodeDecodeError as exc:
            raise ValueError(f'Non-text content in source allowlist: {relative}') from exc
        # Names/public links in the license remain intact; private paths and
        # credentials have no reason to occur there and are checked as elsewhere.
        if (PRIVATE_PATH.search(text) or DELIVERY_ID.search(text) or CREDENTIAL.search(text)
                or any(name in text for name in private_names)):
            raise ValueError(f'Potential private content in {relative}; inspect before exporting')
        snapshots[relative] = (content, stat.S_IMODE(source.stat().st_mode) & 0o777)
    return snapshots


def export(destination, root=None):
    root = Path(root or ROOT).resolve()
    destination = Path(os.path.abspath(os.path.expanduser(str(destination))))
    archive = Path(str(destination) + '.tar.gz')
    for path in (destination, archive):
        if path.exists() or path.is_symlink():
            raise ValueError('Choose a new release output directory and archive; existing contents are never overwritten')
    snapshots = inspect_sources(root)
    manifest = {str(relative): hashlib.sha256(content).hexdigest()
                for relative, (content, mode) in snapshots.items()}
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Build from the inspected snapshots: source edits during packaging cannot
    # introduce unchecked bytes or make the manifest disagree with the archive.
    with tempfile.TemporaryDirectory(prefix='.grip-release-', dir=destination.parent) as temporary:
        stage = Path(temporary) / destination.name
        stage.mkdir()
        for relative, (content, mode) in snapshots.items():
            target = stage / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
            target.chmod(mode)
        (stage / 'RELEASE_MANIFEST.json').write_text(json.dumps({'version': '2.1.0', 'files': manifest}, indent=2) + '\n')
        staged_archive = Path(temporary) / 'source.tar.gz'
        def public_tar_metadata(member):
            member.uid = member.gid = 0
            member.uname = member.gname = ''
            member.mtime = 0
            member.pax_headers = {}
            return member

        # Tar defaults expose the local account/group and timestamps even when
        # file contents are clean. Strip them and omit gzip filename/time fields.
        with staged_archive.open('wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
                with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as tar:
                    tar.add(stage, arcname=destination.name, filter=public_tar_metadata)
        # Exclusive archive creation protects previous releases even if another
        # process creates an archive after the initial existence check.
        with archive.open('xb') as output, staged_archive.open('rb') as source:
            try:
                shutil.copyfileobj(source, output)
                if destination.exists() or destination.is_symlink():
                    raise ValueError('Release output appeared while packaging; refusing to replace it')
                stage.rename(destination)
            except BaseException:
                archive.unlink()
                raise
    result = {'directory': str(destination), 'archive': str(archive), 'files': len(manifest)}
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    export(parser.parse_args().output)
