# workflow/scripts/make_tracks.py
import os

sample = snakemake.wildcards.sample
out_txt = snakemake.output.txt

bw_paths = [
    ("BAM CPM (no blacklist)", snakemake.input.bw1),
    ("BAM CPM (blacklist)", snakemake.input.bw2),
    ("R2 first-base CPM (no blacklist)", snakemake.input.bw3),
    ("R2 first-base CPM (blacklist)", snakemake.input.bw4),
]

os.makedirs(os.path.dirname(out_txt), exist_ok=True)

with open(out_txt, "w") as f:
    for name, path in bw_paths:
        # UCSC track line (works if files are hosted via http(s); for local IGV, you just open bw directly)
        f.write(f'track type=bigWig name="{sample} {name}" description="{sample} {name}" visibility=full autoScale=on bigDataUrl={path}\n')
