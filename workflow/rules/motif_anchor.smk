# Optional motif-conditioned diagnostic subset. Never used as discovery input.
OUTDIR = config["output"]["dir"]
MOTIF_CFG = config.get("motif_anchoring", {})
MOTIF_ENABLED = bool(MOTIF_CFG.get("enable", False))
MOTIF_TARGETS = []
if MOTIF_ENABLED:
    MOTIF_TARGETS = expand(
        f"{OUTDIR}/motif/{{sample}}/{{sample}}.anchored.bed", sample=SAMPLES,
    )

    rule motif_anchor:
        input:
            bed=f"{OUTDIR}/signal/{{sample}}/{{sample}}.crosslink.noblacklist.counts.bed.gz",
            fasta=config["reference"]["fasta"],
            fai=rules.faidx_reference.output.fai
        output:
            bed=f"{OUTDIR}/motif/{{sample}}/{{sample}}.anchored.bed",
            stats=f"{OUTDIR}/motif/{{sample}}/{{sample}}.anchor.stats.json"
        log:
            f"{OUTDIR}/logs/motif/{{sample}}.anchor.log"
        conda:
            "envs/py_signal.yaml"
        params:
            motif=str(MOTIF_CFG.get("motif", "DRACH")),
            motif_pos=int(MOTIF_CFG.get("motif_pos", 3)),
            site_shift=0,  # Input is already an inferred cross-link on RNA strand.
            search_window=int(MOTIF_CFG.get("search_window", 0)),
            keep_unmatched=bool(MOTIF_CFG.get("keep_unmatched", False))
        script:
            "scripts/motif_anchor.py"
