# Independent post-fastp adapter pass: handle residual adapters missed by overlap trimming.
CLEANUP_CFG = config.get("adapter_cleanup", {})
rule adapter_cleanup:
    input:
        r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"
    output:
        r1=temp(f"{OUTDIR}/tmp/adapter_clean/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{OUTDIR}/tmp/adapter_clean/{{sample}}_R2.fastq.gz"),
        json=f"{OUTDIR}/qc/cutadapt/{{sample}}.json",
        audit=f"{OUTDIR}/qc/cutadapt/{{sample}}.residual_audit.json"
    log:
        f"{OUTDIR}/logs/cutadapt/{{sample}}.log"
    threads: 4
    resources:
        mem_mb=3000
    conda:
        "envs/cutadapt.yaml"
    params:
        adapter_r1=CLEANUP_CFG["adapter_r1"],
        adapter_r2=CLEANUP_CFG["adapter_r2"],
        common_prefix=CLEANUP_CFG["exact_common_prefix"] + ";min_overlap=" + str(len(CLEANUP_CFG["exact_common_prefix"])) + ";max_error_rate=0",
        audit_prefix=CLEANUP_CFG["exact_common_prefix"],
        rounds=int(CLEANUP_CFG.get("rounds", 2)),
        error_rate=float(CLEANUP_CFG.get("error_rate", 0.1)),
        overlap=int(CLEANUP_CFG.get("minimum_overlap", 8)),
        minimum_length=int(CLEANUP_CFG.get("minimum_length", 20)),
        audit_script=os.path.join(workflow.basedir, "rules", "scripts", "audit_adapter_fastq.py")
    shell:
        r"""
        mkdir -p $(dirname {output.r1:q}) $(dirname {output.json:q}) $(dirname {log:q})
        cutadapt -j {threads} \
          -a {params.adapter_r1:q} -A {params.adapter_r2:q} \
          -a {params.common_prefix:q} -A {params.common_prefix:q} --times {params.rounds} \
          -e {params.error_rate} -O {params.overlap} \
          -m {params.minimum_length}:{params.minimum_length} --pair-filter=any \
          --json {output.json:q} -o {output.r1:q} -p {output.r2:q} \
          {input.r1:q} {input.r2:q} > {log:q} 2>&1
        python {params.audit_script:q} \
          {wildcards.sample:q} {output.r1:q} {output.r2:q} {output.audit:q} \
          --minimum-length {params.minimum_length} --adapter-prefix {params.audit_prefix:q}
        """
