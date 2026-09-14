"""Run bounded RNA-oriented MEME discovery with locus holdout and deduplication."""
import base64
import html
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

if "snakemake" in globals():
    sys.path.insert(0, str(snakemake.params.script_dir))
from motif_utils import select_discovery_sequences


def run(smk):
    outputs = {key: str(value) for key, value in smk.output.items()}
    label, log_path = str(smk.params.label), str(smk.log[0])
    Path(outputs["html"]).parent.mkdir(parents=True, exist_ok=True)
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    fasta_inputs = [str(smk.input.fasta)] if isinstance(smk.input.fasta, (str, Path)) else list(smk.input.fasta)
    records, stats = select_discovery_sequences(fasta_inputs, int(smk.params.max_sequences),
                                                float(smk.params.holdout_fraction), int(smk.params.seed),
                                                int(smk.params.holdout_block_size), int(smk.params.holdout_guard))
    stats.update({"label": label, "reverse_complement_search": False,
                  "sequence_orientation": "RNA 5prime to 3prime", "independence_unit": "genomic crosslink locus and RNA strand"})
    with open(outputs["sequences"], "w") as out:
        for name, seq in records:
            out.write(f">{name}\n{seq}\n")

    def transparent_logo():
        with open(outputs["logo"], "wb") as out:
            out.write(base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR4nGNgYGBgAAAABQABpfZFQAAAAABJRU5ErkJggg=="))

    if len(records) < int(smk.params.min_sequences):
        stats["status"] = "insufficient_discovery_sequences"
        message = f"MEME not run for {label}: {len(records)} discovery loci; at least {smk.params.min_sequences} required."
        Path(outputs["txt"]).write_text(message + "\n")
        Path(outputs["html"]).write_text("<!doctype html><meta charset='utf-8'><title>GRIP-seq motif</title>" +
                                         f"<h1>{html.escape(label)}</h1><p>{html.escape(message)}</p>")
        transparent_logo()
        Path(log_path).write_text(message + "\n")
    else:
        with tempfile.TemporaryDirectory(prefix="meme_", dir=str(Path(outputs["html"]).parent)) as tmp:
            command = ["meme", outputs["sequences"], "-dna", "-mod", "zoops", "-nmotifs", str(smk.params.nmotifs),
                       "-minw", str(smk.params.min_width), "-maxw", str(smk.params.max_width), "-oc", tmp,
                       "-seed", str(smk.params.seed), "-p", str(smk.threads), "-nostatus"]
            # No -revcomp: the FASTA is already oriented as the original RNA.
            # MEME's width and sequence cap are explicit for reproducible cost.
            with open(log_path, "w") as log:
                log.write("Command: " + " ".join(command) + "\n")
                log.flush()
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            for key, filename in (("html", "meme.html"), ("txt", "meme.txt")):
                shutil.copyfile(os.path.join(tmp, filename), outputs[key])
            logo = os.path.join(tmp, "logo1.png")
            if os.path.exists(logo):
                shutil.copyfile(logo, outputs["logo"])
            else:
                transparent_logo()
            stats["status"] = "completed"
    Path(outputs["stats"]).write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    run(snakemake)
