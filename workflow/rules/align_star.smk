# workflow/rules/align_star.smk
import os

OUTDIR = config["output"]["dir"]

rule star_align_unique:
    input:
        r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"
    output:
        bam=temp(f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"),
        log_final=f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
        log_final_qc=f"{OUTDIR}/qc/star/{{sample}}/{{sample}}.Log.final.out",
        sj=f"{OUTDIR}/star/{{sample}}/{{sample}}.SJ.out.tab"
    log:
        f"{OUTDIR}/logs/star/{{sample}}.log"
    threads: config["threads"]["star"]
    resources:
        mem_mb=40000
    conda:
        "envs/star.yaml"
    params:
        index=config["reference"]["star_index"],
        extra=config["star"].get("extra", ""),
        multimap_nmax=int(config.get("star", {}).get("multimap_nmax", 1)),
        min_match_bases=int(config.get("star", {}).get("min_match_bases", 20)),
        min_score_fraction=float(config.get("star", {}).get("min_score_fraction", 0.0)),
        min_match_fraction=float(config.get("star", {}).get("min_match_fraction", 0.0)),
        max_mismatch_fraction=float(config.get("star", {}).get("max_mismatch_fraction", 0.10)),
        max_mismatches=int(config.get("star", {}).get("max_mismatches", 10))
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {output.log_final}) $(dirname {log})
        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index:q} \
          --readFilesIn {input.r1:q} {input.r2:q} \
          --readFilesCommand zcat \
          --outFilterMultimapNmax {params.multimap_nmax} \
          --outFilterMatchNmin {params.min_match_bases} \
          --outFilterScoreMinOverLread {params.min_score_fraction} \
          --outFilterMatchNminOverLread {params.min_match_fraction} \
          --outFilterMismatchNoverLmax {params.max_mismatch_fraction} \
          --outFilterMismatchNmax {params.max_mismatches} \
          --outFileNamePrefix {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}. \
          --outSAMtype BAM SortedByCoordinate \
          {params.extra} \
          > {log} 2>&1

        # STAR writes:
        # {OUTDIR}/tmp/star/sample/sample.Aligned.sortedByCoord.out.bam
        # {OUTDIR}/tmp/star/sample/sample.Log.final.out
        # {OUTDIR}/tmp/star/sample/sample.SJ.out.tab
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.Log.final.out {output.log_final}
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.Log.final.out {output.log_final_qc}
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.SJ.out.tab {output.sj}
        """


rule summarize_star_mapping:
    input:
        expand(f"{OUTDIR}/qc/star/{{sample}}/{{sample}}.Log.final.out", sample=FASTQ_SAMPLES)
    output:
        f"{OUTDIR}/qc/star/short_read_mapping_summary.tsv"
    script:
        "scripts/summarize_star_logs.py"
