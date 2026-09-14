# workflow/rules/bam_post.smk
import os

OUTDIR = config["output"]["dir"]
MIN_MAPQ = int(config["filtering"]["min_mapq"])

rule filter_mapq_and_index:
    input:
        bam=lambda wc: config["samples"][wc.sample].get("bam", f"{OUTDIR}/tmp/star/{wc.sample}/{wc.sample}.Aligned.sortedByCoord.out.bam")
    output:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai"
    log:
        f"{OUTDIR}/logs/samtools/{{sample}}.mapq_filter.log"
    threads: config["threads"]["samtools"]
    resources:
        mem_mb=4000
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        # Exclude unmapped, secondary, QC-failed and supplementary alignments.
        # PCR/UMI duplicate handling remains an explicit upstream policy.
        samtools quickcheck {input.bam:q}
        samtools view -@ {threads} -b -F 0xB04 -q {MIN_MAPQ} {input.bam:q} \
          | samtools sort -@ {threads} -m 256M -o {output.bam:q} - 2> {log:q}
        samtools index -@ {threads} {output.bam:q} 2>> {log:q}
        """
