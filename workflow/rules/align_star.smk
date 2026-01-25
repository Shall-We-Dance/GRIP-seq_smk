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
        sj=f"{OUTDIR}/star/{{sample}}/{{sample}}.SJ.out.tab"
    log:
        f"logs/star/{{sample}}.log"
    threads: config["threads"]["star"]
    conda:
        "envs/star.yaml"
    params:
        index=config["reference"]["star_index"],
        extra=config["star"].get("extra", "")
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {output.log_final}) $(dirname {log})
        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}. \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattributes NH HI AS nM MD \
          --outFilterType BySJout \
          --outFilterMultimapNmax 1 \
          --outSAMmultNmax 1 \
          --winAnchorMultimapNmax 1 \
          --outFilterMismatchNoverLmax 0.04 \
          --alignSJDBoverhangMin 1 \
          --alignIntronMin 20 \
          --alignIntronMax 1000000 \
          --alignMatesGapMax 1000000 \
          {params.extra} \
          > {log} 2>&1

        # STAR writes:
        # {OUTDIR}/tmp/star/sample/sample.Aligned.sortedByCoord.out.bam
        # {OUTDIR}/tmp/star/sample/sample.Log.final.out
        # {OUTDIR}/tmp/star/sample/sample.SJ.out.tab
        mv {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.Aligned.sortedByCoord.out.bam {output.bam}
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.Log.final.out {output.log_final}
        cp {OUTDIR}/tmp/star/{wildcards.sample}/{wildcards.sample}.SJ.out.tab {output.sj}
        """
