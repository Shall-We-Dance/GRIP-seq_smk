"""Average replicate CPM tracks with missing bases equal to zero."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(snakemake.scriptdir) if "snakemake" in globals() else str(Path(__file__).resolve().parent))
from bigwig_math import combine_tracks


def main():
    if "snakemake" in globals():
        inputs, output = list(snakemake.input.bws), snakemake.output.bw
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("--bws", nargs="+", required=True)
        parser.add_argument("--out-bw", required=True)
        args = parser.parse_args()
        inputs, output = args.bws, args.out_bw
    combine_tracks(inputs, output)


if __name__ == "__main__":
    main()
