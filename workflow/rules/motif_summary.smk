# Include after motif_de_novo.smk and motif_distribution.smk. Summaries consume
# only completed per-sample results; no sites are pooled for new tests.
MOTIF_SUMMARY_TARGETS = []
if MOTIF_DIST_ENABLED:
    MOTIF_SUMMARY_TARGETS = [
        f"{OUTDIR}/motif_distribution/integrated.summary.tsv",
        f"{OUTDIR}/motif_distribution/grouped.predefined_motifs.pdf",
        f"{OUTDIR}/motif_distribution/grouped.curves.tsv",
    ]

    rule summarize_motif_offsets:
        input:
            tables=expand(
                f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.offsets.tsv",
                sample=MOTIF_DIST_SAMPLES,
            ),
            statistics=expand(
                f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.statistics.json",
                sample=MOTIF_DIST_SAMPLES,
            )
        output:
            summary=f"{OUTDIR}/motif_distribution/integrated.summary.tsv",
            plot=f"{OUTDIR}/motif_distribution/grouped.predefined_motifs.pdf",
            curves=f"{OUTDIR}/motif_distribution/grouped.curves.tsv"
        params:
            samples=MOTIF_DIST_SAMPLES,
            groups=config.get("bigwigCompare", {}).get("groups", {})
        resources:
            mem_mb=1500
        conda:
            "envs/qc.yaml"
        script:
            "scripts/summarize_motif_offsets.py"
