# workflow/scripts/make_tracks.py
import os

sample = snakemake.wildcards.sample
out_txt = snakemake.output.txt

def add_track(tracks, name, path):
    if not path:
        return
    if isinstance(path, (list, tuple)):
        if not path:
            return
        path = path[0]
    tracks.append((name, str(path)))


bw_paths = []
add_track(bw_paths, "BAM CPM (no blacklist)", snakemake.input.get("bw1"))
add_track(bw_paths, "BAM CPM (blacklist)", snakemake.input.get("bw2"))
add_track(bw_paths, "R2 first-base CPM (no blacklist)", snakemake.input.get("bw3"))
add_track(bw_paths, "R2 first-base CPM (blacklist)", snakemake.input.get("bw4"))

os.makedirs(os.path.dirname(out_txt), exist_ok=True)

with open(out_txt, "w") as f:
    for name, path in bw_paths:
        # UCSC track line (works if files are hosted via http(s); for local IGV, you just open bw directly)
        f.write(f'track type=bigWig name="{sample} {name}" description="{sample} {name}" visibility=full autoScale=on bigDataUrl={path}\n')
