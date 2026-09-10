OUTDIR = config["output"]["dir"]
METAPLOT_CFG = config.get("metaplot", {})
METAPLOT_ENABLED = bool(METAPLOT_CFG.get("enable", False))
configured_metaplot_samples = METAPLOT_CFG.get("samples", "all")
METAPLOT_SAMPLES = SAMPLES if configured_metaplot_samples == "all" else list(configured_metaplot_samples)
METAPLOT_FORMAT = str(METAPLOT_CFG.get("plot_format", "pdf"))
METAPLOT_TRACKS = list(METAPLOT_CFG.get(
    "tracks", ["bamCPM.noblacklist", "R2firstbaseCPM.noblacklist"]
))

METAPLOT_TARGETS = []
if METAPLOT_ENABLED:
    METAPLOT_TARGETS = expand(
        f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.TSS_TES.{METAPLOT_FORMAT}",
        sample=METAPLOT_SAMPLES,
    ) + [f"{OUTDIR}/metaplot/integrated.TSS_TES.{METAPLOT_FORMAT}"]


def metaplot_bigwigs(wildcards):
    return [
        f"{OUTDIR}/bigwig/{wildcards.sample}/{wildcards.sample}.{track}.bw"
        for track in METAPLOT_TRACKS
    ]


if METAPLOT_ENABLED:
    rule prepare_metagene_regions:
        input:
            gtf=config["reference"]["gtf"]
        output:
            bed=f"{OUTDIR}/metaplot/reference.genes.bed"
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/gtf_to_gene_bed.py"

    rule compute_metaplot_matrix:
        input:
            bws=metaplot_bigwigs,
            bed=rules.prepare_metagene_regions.output.bed
        output:
            matrix=temp(f"{OUTDIR}/tmp/metaplot/{{sample}}.matrix.gz")
        log:
            "logs/metaplot/{sample}.computeMatrix.log"
        threads:
            int(config.get("threads", {}).get("bigwig", 8))
        conda:
            "envs/qc.yaml"
        params:
            upstream=int(METAPLOT_CFG.get("upstream", 3000)),
            downstream=int(METAPLOT_CFG.get("downstream", 3000)),
            body=int(METAPLOT_CFG.get("gene_body_length", 5000)),
            bin_size=int(METAPLOT_CFG.get("bin_size", 50))
        shell:
            r"""
            mkdir -p $(dirname {output.matrix}) $(dirname {log})
            computeMatrix scale-regions \
              -S {input.bws} -R {input.bed} \
              --beforeRegionStartLength {params.upstream} \
              --regionBodyLength {params.body} \
              --afterRegionStartLength {params.downstream} \
              --binSize {params.bin_size} \
              --averageTypeBins mean --missingDataAsZero \
              --numberOfProcessors {threads} \
              -o {output.matrix} > {log} 2>&1
            """

    rule plot_metaplot_sample:
        input:
            matrix=f"{OUTDIR}/tmp/metaplot/{{sample}}.matrix.gz"
        output:
            plot=f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.TSS_TES.{METAPLOT_FORMAT}",
            data=f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.TSS_TES.profile.tsv"
        log:
            "logs/metaplot/{sample}.plotProfile.log"
        conda:
            "envs/qc.yaml"
        params:
            labels=" ".join(f"'{track}'" for track in METAPLOT_TRACKS),
            title=lambda wc: f"'{wc.sample}: TSS-TES (strand aware)'"
        shell:
            r"""
            mkdir -p $(dirname {output.plot}) $(dirname {log})
            plotProfile -m {input.matrix} -out {output.plot} \
              --samplesLabel {params.labels} \
              --plotTitle {params.title} \
              --regionsLabel genes \
              --plotHeight 6 --plotWidth 9 \
              --outFileNameData {output.data} > {log} 2>&1
            """

    rule plot_metaplot_integrated:
        input:
            expand(
                f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.TSS_TES.profile.tsv",
                sample=METAPLOT_SAMPLES,
            )
        output:
            plot=f"{OUTDIR}/metaplot/integrated.TSS_TES.{METAPLOT_FORMAT}"
        conda:
            "envs/qc.yaml"
        params:
            samples=METAPLOT_SAMPLES,
            upstream=int(METAPLOT_CFG.get("upstream", 3000)),
            downstream=int(METAPLOT_CFG.get("downstream", 3000)),
            body=int(METAPLOT_CFG.get("gene_body_length", 5000)),
            bin_size=int(METAPLOT_CFG.get("bin_size", 50)),
            n_bins=(
                int(METAPLOT_CFG.get("upstream", 3000))
                + int(METAPLOT_CFG.get("gene_body_length", 5000))
                + int(METAPLOT_CFG.get("downstream", 3000))
            ) // int(METAPLOT_CFG.get("bin_size", 50))
        script:
            "scripts/plot_combined_metaplot.py"
