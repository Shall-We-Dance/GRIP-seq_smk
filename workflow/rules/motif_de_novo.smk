OUTDIR = config["output"]["dir"]
DE_NOVO_CFG = config.get("de_novo_motif", {})
DE_NOVO_ENABLED = bool(DE_NOVO_CFG.get("enable", False))
configured_motif_samples = DE_NOVO_CFG.get("samples", "all")
MOTIF_SAMPLES = SAMPLES if configured_motif_samples == "all" else list(configured_motif_samples)

DE_NOVO_MOTIF_TARGETS = []
if DE_NOVO_ENABLED:
    DE_NOVO_MOTIF_TARGETS = expand(
        f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.meme.html",
        sample=MOTIF_SAMPLES,
    ) + expand(
        f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.motif1.png",
        sample=MOTIF_SAMPLES,
    ) + [
        f"{OUTDIR}/motif_de_novo/integrated/integrated.meme.html",
        f"{OUTDIR}/motif_de_novo/integrated/integrated.motif1.png",
    ]

if DE_NOVO_ENABLED:
    rule call_rt_termination_sites:
        input:
            bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
            bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
            fasta=config["reference"]["fasta"],
            fai=rules.faidx_reference.output.fai
        output:
            bed=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.bed",
            fasta=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.fa",
            stats=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.stats.tsv"
        log:
            "logs/motif_de_novo/{sample}.sites.log"
        conda:
            "envs/py_signal.yaml"
        params:
            min_mapq=int(config["filtering"]["min_mapq"]),
            min_reads=int(DE_NOVO_CFG.get("min_reads", 40)),
            neighbor_fold=float(DE_NOVO_CFG.get("neighbor_fold", 1.5)),
            flank=int(DE_NOVO_CFG.get("sequence_flank", 10))
        script:
            "scripts/call_rt_termination_sites.py"

    rule meme_de_novo_sample:
        input:
            fasta=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.fa"
        output:
            html=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.meme.html",
            txt=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.meme.txt",
            logo=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.motif1.png"
        log:
            "logs/motif_de_novo/{sample}.meme.log"
        conda:
            "envs/motif.yaml"
        params:
            label=lambda wc: wc.sample,
            min_sequences=int(DE_NOVO_CFG.get("min_sequences", 10)),
            nmotifs=int(DE_NOVO_CFG.get("nmotifs", 5)),
            min_width=int(DE_NOVO_CFG.get("min_width", 4)),
            max_width=int(DE_NOVO_CFG.get("max_width", 8))
        script:
            "scripts/run_meme.py"

    rule combine_motif_sequences:
        input:
            expand(
                f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.fa",
                sample=MOTIF_SAMPLES,
            )
        output:
            fasta=f"{OUTDIR}/motif_de_novo/integrated/integrated.rt_termination.fa"
        shell:
            "mkdir -p $(dirname {output.fasta}) && cat {input} > {output.fasta}"

    rule meme_de_novo_integrated:
        input:
            fasta=rules.combine_motif_sequences.output.fasta
        output:
            html=f"{OUTDIR}/motif_de_novo/integrated/integrated.meme.html",
            txt=f"{OUTDIR}/motif_de_novo/integrated/integrated.meme.txt",
            logo=f"{OUTDIR}/motif_de_novo/integrated/integrated.motif1.png"
        log:
            "logs/motif_de_novo/integrated.meme.log"
        conda:
            "envs/motif.yaml"
        params:
            label="integrated IP samples",
            min_sequences=int(DE_NOVO_CFG.get("min_sequences", 10)),
            nmotifs=int(DE_NOVO_CFG.get("nmotifs", 5)),
            min_width=int(DE_NOVO_CFG.get("min_width", 4)),
            max_width=int(DE_NOVO_CFG.get("max_width", 8))
        script:
            "scripts/run_meme.py"
