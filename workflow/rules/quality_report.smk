# Runtime and library QC stays local; packaged releases omit all run artifacts.
QC_TARGETS = [f"{OUTDIR}/qc/library_qc.tsv", f"{OUTDIR}/qc/library_qc.json", f"{OUTDIR}/qc/library_qc.pdf", f"{OUTDIR}/provenance.json"]

rule library_quality_report:
    input:
        signals=expand(f"{OUTDIR}/signal/{{sample}}/{{sample}}.noblacklist.signal.stats.json", sample=SAMPLES),
        gene_signals=expand(f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.gene_signal.tsv", sample=METAPLOT_SAMPLES) if METAPLOT_ENABLED else []
    output:
        table=f"{OUTDIR}/qc/library_qc.tsv",
        json=f"{OUTDIR}/qc/library_qc.json",
        plot=f"{OUTDIR}/qc/library_qc.pdf",
        correlations=f"{OUTDIR}/qc/gene_signal_correlations.tsv"
    conda:
        "envs/qc.yaml"
    script:
        "scripts/library_quality_report.py"

rule analysis_provenance:
    input:
        qc=rules.library_quality_report.output.table,
        bams=expand(f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam", sample=SAMPLES)
    output:
        f"{OUTDIR}/provenance.json"
    script:
        "scripts/write_provenance.py"
