# workflow/rules/bam_post.smk
import os

OUTDIR = config["output"]["dir"]
MIN_MAPQ = int(config["filtering"]["min_mapq"])

rule filter_mapq_and_index:
    input:
        bam=f"{OUTDIR}/tmp/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai"
    log:
        f"logs/samtools/{{sample}}.mapq_filter.log"
    threads: config["threads"]["samtools"]
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        # MAPQ > 10  => min_mapq=11; drop secondary (0x100) and supplementary (0x800)
        samtools view -@ {threads} -b -F 0x100 -F 0x800 -q {MIN_MAPQ} {input.bam} \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """
