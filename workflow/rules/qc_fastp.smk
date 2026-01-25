# workflow/rules/qc_fastp.smk
import os

OUTDIR = config["output"]["dir"]

def units(sample):
    return list(range(len(config["samples"][sample]["R1"])))

rule fastp_unit:
    input:
        r1=lambda wc: config["samples"][wc.sample]["R1"][int(wc.unit)],
        r2=lambda wc: config["samples"][wc.sample]["R2"][int(wc.unit)]
    output:
        clean_r1=temp(f"{OUTDIR}/tmp/fastp/{wildcards.sample}/unit{wildcards.unit}_R1.fastq.gz"),
        clean_r2=temp(f"{OUTDIR}/tmp/fastp/{wildcards.sample}/unit{wildcards.unit}_R2.fastq.gz"),
        html=f"{OUTDIR}/qc/fastp/{{sample}}/unit{{unit}}.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}/unit{{unit}}.json"
    log:
        f"logs/fastp/{{sample}}/unit{{unit}}.log"
    threads: config["threads"]["fastp"]
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.clean_r1}) $(dirname {output.html}) $(dirname {log})
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.clean_r1} -O {output.clean_r2} \
          --html {output.html} --json {output.json} \
          --thread {threads} \
          > {log} 2>&1
        """

rule merge_fastq_after_fastp:
    input:
        r1=lambda wc: [f"{OUTDIR}/tmp/fastp/{wc.sample}/unit{i}_R1.fastq.gz" for i in units(wc.sample)],
        r2=lambda wc: [f"{OUTDIR}/tmp/fastp/{wc.sample}/unit{i}_R2.fastq.gz" for i in units(wc.sample)],
        html=lambda wc: [f"{OUTDIR}/qc/fastp/{wc.sample}/unit{i}.html" for i in units(wc.sample)],
        json=lambda wc: [f"{OUTDIR}/qc/fastp/{wc.sample}/unit{i}.json" for i in units(wc.sample)],
    output:
        merged_r1=temp(f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R1.fastq.gz"),
        merged_r2=temp(f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R2.fastq.gz")
    threads: 2
    shell:
        r"""
        mkdir -p $(dirname {output.merged_r1})
        cat {input.r1} > {output.merged_r1}
        cat {input.r2} > {output.merged_r2}
        """

rule fastp_sample_report_only:
    input:
        r1=f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R2.fastq.gz"
    output:
        # sample-level report (kept)
        html=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.json",
        # sample-level output fastq (temp, not kept)
        out_r1=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz"),
        out_r2=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"),
    log:
        f"logs/fastp/{{sample}}/sample_report_only.log"
    threads: max(2, int(config["threads"]["fastp"]))
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.html}) $(dirname {output.out_r1}) $(dirname {log})
        # "report-only" in practice: disable trimming/filtering so reads are not changed materially
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.out_r1} -O {output.out_r2} \
          --disable_adapter_trimming \
          --disable_quality_filtering \
          --disable_length_filtering \
          --html {output.html} --json {output.json} \
          --thread {threads} \
          > {log} 2>&1
        """

rule multiqc:
    input:
        # depend on sample-level fastp report for each sample (unit reports will be picked up too)
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html", sample=list(config["samples"].keys()))
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc -o {OUTDIR}/qc/multiqc {OUTDIR} > {log} 2>&1
        """
