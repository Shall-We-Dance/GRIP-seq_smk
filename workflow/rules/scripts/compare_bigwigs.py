"""Exact 1-bp comparison: same per-base arithmetic as bigwigCompare."""
from pathlib import Path
import sys
sys.path.insert(0, str(snakemake.scriptdir) if "snakemake" in globals() else str(Path(__file__).resolve().parent))
from bigwig_math import combine_tracks
combine_tracks([snakemake.input.bw1, snakemake.input.bw2], snakemake.output.bw,
               snakemake.params.operation, float(snakemake.params.pseudocount))
