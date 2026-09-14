"""Prepare a mature-transcript cohort for segmented UTR/CDS metaplots.

The historical filename is retained for workflow compatibility; outputs are
JSONL block annotations plus explicit selection/QC tables, not gene-span BED.
"""

import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(snakemake.scriptdir))
from metagene import common_chromosome_sizes, prepare_annotation
import pyBigWig

headers = []
for path in dict.fromkeys(snakemake.input.bws):
    with pyBigWig.open(str(path)) as bigwig:
        headers.append(bigwig.chroms())
common_sizes, reference_summary = common_chromosome_sizes(headers)
models, rows, summary = prepare_annotation(
    snakemake.input.gtf, int(snakemake.params.min_region_length), common_sizes
)
summary["bigwig_reference_cohort"] = reference_summary
for path in snakemake.output:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.models, "w") as handle:
    for model in models:
        handle.write(json.dumps(model, sort_keys=True) + "\n")
with open(snakemake.output.selection, "w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
with open(snakemake.output.summary, "w") as handle:
    json.dump(summary, handle, indent=2, sort_keys=True)
    handle.write("\n")
