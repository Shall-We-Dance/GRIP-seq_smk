"""Write local track descriptions without uploading any data."""
from pathlib import Path
sample = str(snakemake.wildcards.sample)
output = Path(snakemake.output.txt)
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as f:
    for path in snakemake.input.tracks:
        label = Path(path).name.removeprefix(sample + ".").removesuffix(".bw")
        f.write(f'track type=bigWig name="{sample} {label}" description="{sample} {label}" visibility=full autoScale=on bigDataUrl={path}\n')
