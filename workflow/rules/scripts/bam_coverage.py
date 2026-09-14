"""Create aligned-read CPM coverage or copy an externally validated baseline.

Reuse is an explicit BAM-mode option. The caller must establish matching BAM
selection, reference, blacklist and bamCoverage parameters before configuring it.
Outputs are independent file copies, so later reruns cannot alter source tracks.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def create_coverage(sm):
    destination = Path(sm.output.bw)
    destination.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(sm.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    reuse = list(sm.input.get('reuse', []))
    if len(reuse) > 1:
        raise ValueError('At most one validated coverage source is allowed')
    with log_path.open('w') as log:
        if reuse:
            source = Path(reuse[0])
            if source.resolve() == destination.resolve():
                raise ValueError('Coverage reuse source must differ from its output')
            import pysam
            import pyBigWig
            with pysam.AlignmentFile(str(sm.input.bam), 'rb') as bam:
                expected_chroms = dict(zip(bam.references, bam.lengths))
            with pyBigWig.open(str(source)) as track:
                if not track.isBigWig() or track.chroms() != expected_chroms:
                    raise ValueError('Reused coverage BigWig must match BAM chromosome names and lengths')
            # Stage beside the output and replace its directory entry atomically.
            # This also breaks pre-existing links instead of overwriting targets.
            descriptor, staged = tempfile.mkstemp(prefix='.coverage-', suffix='.bw', dir=destination.parent)
            os.close(descriptor)
            try:
                shutil.copyfile(source, staged)
                os.chmod(staged, source.stat().st_mode & 0o777)
                os.replace(staged, destination)
            finally:
                if os.path.exists(staged):
                    os.unlink(staged)
            json.dump({'mode': 'reused_validated_bam_coverage', 'source': str(source),
                       'bam_dependency': str(sm.input.bam),
                       'blacklist_dependency': str(sm.input.get('bl', '')),
                       'copy_policy': 'independent bytes; no symlink or hardlink'}, log, indent=2)
            log.write('\n')
            return
        command = ['bamCoverage', '-b', str(sm.input.bam), '-o', str(destination),
                   '--binSize', '1', '--normalizeUsing', 'CPM',
                   '--numberOfProcessors', str(sm.threads)]
        if sm.input.get('bl'):
            command += ['--blackListFileName', str(sm.input.bl)]
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


if __name__ == '__main__':
    create_coverage(snakemake)
