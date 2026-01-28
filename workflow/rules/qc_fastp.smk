# workflow/rules/qc_fastp.smk
import os

OUTDIR = config["output"]["dir"]

def units(sample):
    return list(range(len(config["samples"][sample]["R1"])))


# ----------------------------
# Per-lane (unit) fastp
#   Step1: QC + Illumina adapter trimming + dedup
#   Step2: GRIP-seq custom trimming (fixed trims; disable adapter trimming)
# Outputs:
#   - cleaned FASTQ (temp) for merge
#   - final per-lane fastp HTML/JSON (kept; corresponds to trimmed reads used downstream)
# ----------------------------
rule fastp_unit:
    input:
        r1=lambda wc: config["samples"][wc.sample]["R1"][int(wc.unit)],
        r2=lambda wc: config["samples"][wc.sample]["R2"][int(wc.unit)]
    output:
        clean_r1=temp(f"{OUTDIR}/tmp/fastp/{{sample}}/unit{{unit}}_R1.fastq.gz"),
        clean_r2=temp(f"{OUTDIR}/tmp/fastp/{{sample}}/unit{{unit}}_R2.fastq.gz"),
        html_step1=f"{OUTDIR}/qc/fastp/{{sample}}/unit{{unit}}_step1.html",
        json_step1=f"{OUTDIR}/qc/fastp/{{sample}}/unit{{unit}}_step1.json",
        html_step2=f"{OUTDIR}/qc/fastp/{{sample}}/unit{{unit}}_step2.html",
        json_step2=f"{OUTDIR}/qc/fastp/{{sample}}/unit{{unit}}_step2.json"
    log:
        f"logs/fastp/{{sample}}/unit{{unit}}.log"
    threads: int(config["threads"]["fastp"])
    conda:
        "envs/qc.yaml"
    params:
        dedup_arg=("--dedup" if bool(config.get("fastp", {}).get("dedup_adapter", {}).get("dedup", True)) else ""),
        disable_adapter_arg=("--disable_adapter_trimming" if bool(config.get("fastp", {}).get("grip_trim", {}).get("disable_adapter_trimming", True)) else ""),
        trim_f=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_front_r1", 12)),
        trim_t=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_tail_r1", 10)),
        trim_F=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_front_r2", 10)),
        trim_T=int(config.get("fastp", {}).get("grip_trim", {}).get("trim_tail_r2", 12)),
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.clean_r1}) $(dirname {output.html_step1}) $(dirname {log})

        TMP_R1={OUTDIR}/tmp/fastp/{wildcards.sample}/unit{wildcards.unit}.step1_R1.fastq.gz
        TMP_R2={OUTDIR}/tmp/fastp/{wildcards.sample}/unit{wildcards.unit}.step1_R2.fastq.gz

        # Step 1: QC + Illumina adapter trimming + dedup
        fastp \
          -i {input.r1} -I {input.r2} \
          -o $TMP_R1 -O $TMP_R2 \
          --thread {threads} \
          {params.dedup_arg} \
          --html {output.html_step1} --json {output.json_step1} \
          >> {log} 2>&1

        # Step 2: GRIP-seq custom trimming (disable adapter trimming + fixed trims)
        fastp \
          -i $TMP_R1 -I $TMP_R2 \
          -o {output.clean_r1} -O {output.clean_r2} \
          --thread {threads} \
          {params.disable_adapter_arg} \
          -f {params.trim_f} -t {params.trim_t} \
          -F {params.trim_F} -T {params.trim_T} \
          --html {output.html_step2} --json {output.json_step2} \
          >> {log} 2>&1
        """


# ----------------------------
# Merge post-fastp per sample (after per-lane processing)
# ----------------------------
rule merge_fastq_after_fastp:
    input:
        r1=lambda wc: [f"{OUTDIR}/tmp/fastp/{wc.sample}/unit{i}_R1.fastq.gz" for i in units(wc.sample)],
        r2=lambda wc: [f"{OUTDIR}/tmp/fastp/{wc.sample}/unit{i}_R2.fastq.gz" for i in units(wc.sample)],
        # Ensure per-lane reports exist before merging (optional but helpful for DAG clarity)
        html=lambda wc: [f"{OUTDIR}/qc/fastp/{wc.sample}/unit{i}_step2.html" for i in units(wc.sample)],
        json=lambda wc: [f"{OUTDIR}/qc/fastp/{wc.sample}/unit{i}_step2.json" for i in units(wc.sample)],
    output:
        merged_r1=temp(f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R1.fastq.gz"),
        merged_r2=temp(f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R2.fastq.gz")
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged_r1})
        cat {input.r1} > {output.merged_r1}
        cat {input.r2} > {output.merged_r2}
        """


# ----------------------------
# Sample-level fastp report-only (do not modify reads)
#   Uses merged FASTQs, generates sample-level HTML/JSON.
#   Output FASTQs are temp and not kept.
# ----------------------------
rule fastp_sample_report_only:
    input:
        r1=f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/merged_fastq/{{sample}}_R2.fastq.gz"
    output:
        html=f"{OUTDIR}/qc/fastp/{{sample}}/merged_fastp_final.html",
        json=f"{OUTDIR}/qc/fastp/{{sample}}/merged_fastp_final.json",
        out_r1=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz"),
        out_r2=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"),
    log:
        f"logs/fastp/{{sample}}/sample_report_only.log"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {output.out_r1}) $(dirname {log})

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


# ----------------------------
# MultiQC summary
#   Use --force to avoid _1 suffix when outputs already exist.
#   Limit scanning scope to QC + STAR outputs for performance/stability.
# ----------------------------
rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/merged_fastp_final.html", sample=list(config["samples"].keys()))
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc --force -o {OUTDIR}/qc/multiqc {OUTDIR}/qc {OUTDIR}/star > {log} 2>&1
        """
