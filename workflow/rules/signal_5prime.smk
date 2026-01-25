# workflow/rules/signal_5prime.smk  (part 1: BAM -> BigWig)
import os

OUTDIR = config["output"]["dir"]

rule bam_to_bigwig_cpm_noblacklist:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai"
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.noblacklist.bw"
    log:
        f"logs/bigwig/{{sample}}.bamCoverage.noblacklist.log"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --binSize 1 \
          --normalizeUsing CPM \
          --numberOfProcessors {threads} \
          > {log} 2>&1
        """

rule bam_to_bigwig_cpm_blacklist:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
        bl=rules.get_blacklist.output.bed
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.blacklist.bw"
    log:
        f"logs/bigwig/{{sample}}.bamCoverage.blacklist.log"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        bamCoverage \
          -b {input.bam} \
          -o {output.bw} \
          --binSize 1 \
          --normalizeUsing CPM \
          --blackListFileName {input.bl} \
          --numberOfProcessors {threads} \
          > {log} 2>&1
        """

# workflow/rules/signal_5prime.smk  (part 2: Read2-first-base -> BigWig)
import os

OUTDIR = config["output"]["dir"]

rule r2_first_base_cpm_noblacklist:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
        fai=rules.faidx_reference.output.fai
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.noblacklist.bw"
    log:
        f"logs/5prime/{{sample}}.R2firstbase.noblacklist.log"
    conda:
        "envs/py_signal.yaml"
    params:
        min_mapq=int(config["filtering"]["min_mapq"])
    shell:
        r"""
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        python scripts/extract_r2_first_base.py \
          --bam {input.bam} \
          --fai {input.fai} \
          --out-bw {output.bw} \
          --min-mapq {params.min_mapq} \
          > {log} 2>&1
        """

rule r2_first_base_cpm_blacklist:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
        fai=rules.faidx_reference.output.fai,
        bl=rules.get_blacklist.output.bed
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.blacklist.bw"
    log:
        f"logs/5prime/{{sample}}.R2firstbase.blacklist.log"
    conda:
        "envs/py_signal.yaml"
    params:
        min_mapq=int(config["filtering"]["min_mapq"])
    shell:
        r"""
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        python scripts/extract_r2_first_base.py \
          --bam {input.bam} \
          --fai {input.fai} \
          --out-bw {output.bw} \
          --min-mapq {params.min_mapq} \
          --blacklist-bed {input.bl} \
          > {log} 2>&1
        """
