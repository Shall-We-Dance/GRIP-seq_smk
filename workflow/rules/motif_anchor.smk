# workflow/rules/motif_anchor.smk
import os

OUTDIR = config["output"]["dir"]
MOTIF_CFG = config.get("motif_anchoring", {})
MOTIF_ENABLED = bool(MOTIF_CFG.get("enable", True))


def motif_anchor_targets():
    if not MOTIF_ENABLED:
        return []
    return [
        f"{OUTDIR}/motif/{{sample}}/{{sample}}.anchored.bed",
    ]


MOTIF_TARGETS = [t.format(sample=s) for s in SAMPLES for t in motif_anchor_targets()]


if MOTIF_ENABLED:
    # Strand-aware 1-nt positions of R2 first base (basis for motif anchoring)
    rule r2_first_base_positions_bed:
        input:
            bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
            bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
            fai=rules.faidx_reference.output.fai
        output:
            bed=temp(f"{OUTDIR}/tmp/motif/{{sample}}/{{sample}}.R2firstbase.pos.bed")
        log:
            f"logs/motif/{{sample}}.positions.log"
        conda:
            "envs/py_signal.yaml"
        params:
            min_mapq=int(config["filtering"]["min_mapq"])
        script:
            "scripts/extract_r2_first_base.py"

    # Optional motif anchoring (default DRACH, m6A = 3rd base A); output corrected 1-nt BED
    rule motif_anchor:
        input:
            bed=f"{OUTDIR}/tmp/motif/{{sample}}/{{sample}}.R2firstbase.pos.bed",
            fasta=config["reference"]["fasta"],
            fai=rules.faidx_reference.output.fai
        output:
            bed=f"{OUTDIR}/motif/{{sample}}/{{sample}}.anchored.bed"
        log:
            f"logs/motif/{{sample}}.anchor.log"
        conda:
            "envs/py_signal.yaml"
        params:
            motif=str(MOTIF_CFG.get("motif", "DRACH")),
            motif_pos=int(MOTIF_CFG.get("motif_pos", 3)),
            site_shift=int(MOTIF_CFG.get("site_shift", 1)),
            search_window=int(MOTIF_CFG.get("search_window", 0)),
            keep_unmatched=bool(MOTIF_CFG.get("keep_unmatched", False))
        script:
            "scripts/motif_anchor.py"
