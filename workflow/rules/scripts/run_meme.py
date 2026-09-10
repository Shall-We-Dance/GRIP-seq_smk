"""Run MEME, producing explicit outputs even for low-yield libraries."""

import html
import base64
import os
import shutil
import subprocess
import tempfile


fasta = str(snakemake.input.fasta)
label = str(snakemake.params.label)
min_sequences = int(snakemake.params.min_sequences)
outputs = {key: str(value) for key, value in snakemake.output.items()}
log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(outputs["html"]), exist_ok=True)
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(fasta) as handle:
    n_sequences = sum(1 for line in handle if line.startswith(">"))


def write_transparent_logo():
    """Write a valid placeholder without touching a completed MEME report."""
    transparent_png = (
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVQIHWP4"
        "z8DwHwAFgAI/ScLJAAAAAElFTkSuQmCC"
    )
    with open(outputs["logo"], "wb") as handle:
        handle.write(base64.b64decode(transparent_png))


def placeholder(message):
    with open(outputs["txt"], "w") as handle:
        handle.write(message + "\n")
    with open(outputs["html"], "w") as handle:
        handle.write(
            "<!doctype html><meta charset='utf-8'><title>GRIP-seq motif</title>"
            f"<h1>{html.escape(label)}</h1><p>{html.escape(message)}</p>"
        )
    write_transparent_logo()


if n_sequences < min_sequences:
    message = (
        f"MEME not run for {label}: {n_sequences} sequences; "
        f"at least {min_sequences} required."
    )
    placeholder(message)
    with open(log_path, "w") as log:
        log.write(message + "\n")
else:
    with tempfile.TemporaryDirectory(prefix="meme_", dir=os.path.dirname(outputs["html"])) as tmp:
        command = [
            "meme",
            fasta,
            "-dna",
            "-mod",
            "zoops",
            "-nmotifs",
            str(snakemake.params.nmotifs),
            "-minw",
            str(snakemake.params.min_width),
            "-maxw",
            str(snakemake.params.max_width),
            "-oc",
            tmp,
            "-nostatus",
        ]
        with open(log_path, "w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        shutil.copyfile(os.path.join(tmp, "meme.html"), outputs["html"])
        shutil.copyfile(os.path.join(tmp, "meme.txt"), outputs["txt"])
        logo = os.path.join(tmp, "logo1.png")
        if os.path.exists(logo):
            shutil.copyfile(logo, outputs["logo"])
        else:
            write_transparent_logo()
            with open(log_path, "a") as log:
                log.write(f"MEME completed for {label}, but motif-1 logo was not generated.\n")
