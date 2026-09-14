"""Private runtime provenance; distribution packages exclude this output."""
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


root=Path(snakemake.scriptdir).resolve().parents[2]
inputs=[]
for sample,settings in snakemake.config['samples'].items():
    paths=[settings['bam']] if settings.get('bam') else settings['R1']+settings['R2']
    for path in paths:
        p=Path(path);st=p.stat()
        record={'sample':sample,'path':str(p.resolve()),'size':st.st_size,'mtime_ns':st.st_mtime_ns}
        # BAM reanalysis inputs are compact enough for full checksums.
        if settings.get('bam'): record['sha256']=digest(p)
        inputs.append(record)
coverage_inputs=[]
for sample,settings in snakemake.config['samples'].items():
    for variant,path in settings.get('bam_coverage',{}).items():
        coverage_inputs.append({'sample':sample,'variant':variant,'path':str(Path(path).resolve()),'sha256':digest(path)})
code={str(p.relative_to(root)):digest(p) for p in sorted((root/'workflow').rglob('*')) if p.is_file() and p.suffix in {'.py','.smk','.yaml'} or p.is_file() and p.name=='Snakefile'}
reference_inputs=[]
for kind in ('fasta','gtf'):
    value=snakemake.config.get('reference',{}).get(kind)
    if value and Path(value).is_file():
        reference_inputs.append({'kind':kind,'path':str(Path(value).resolve()),'sha256':digest(value)})
environments={}
for meta in sorted((root/'.snakemake/conda').glob('*_/conda-meta')):
    packages=[]
    for record in sorted(meta.glob('*.json')):
        package=json.loads(record.read_text())
        packages.append({key:package.get(key) for key in ('name','version','build')})
    environments[meta.parent.name]=packages
obj={'schema_version':1,'workflow_version':'2.0.0','created_utc':datetime.now(timezone.utc).isoformat(),
     'python':platform.python_version(),'reference_inputs':reference_inputs,'installed_rule_environments':environments,'config':dict(snakemake.config),'inputs':inputs,'imported_coverage':coverage_inputs,'code_sha256':code,
     'coordinate_system':'0-based half-open; endpoint and inferred crosslink are distinct; all site BED strands denote RNA',
     'reanalysis_note':snakemake.config.get('provenance',{})}
validation=snakemake.config.get('provenance',{}).get('coverage_reuse_validation')
if validation:
    obj['coverage_reuse_validation']=json.loads(Path(validation).read_text())
Path(snakemake.output[0]).write_text(json.dumps(obj,indent=2,default=str)+'\n')
