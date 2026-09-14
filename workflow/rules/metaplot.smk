OUTDIR = config["output"]["dir"]
METAPLOT_CFG = config.get("metaplot", {})
METAPLOT_ENABLED = bool(METAPLOT_CFG.get("enable", False))
configured_metaplot_samples = METAPLOT_CFG.get("samples", "all")
METAPLOT_SAMPLES = SAMPLES if configured_metaplot_samples == "all" else list(configured_metaplot_samples)
METAPLOT_FORMAT = str(METAPLOT_CFG.get("plot_format", "pdf"))
METAPLOT_TRACKS = list(METAPLOT_CFG.get(
    "tracks", ["R2firstbaseCPM.noblacklist", "crosslinkCPM.noblacklist"]
))
METAPLOT_BINS = METAPLOT_CFG.get("bins", {"utr5": 100, "cds": 100, "utr3": 100})
METAPLOT_AXIS_SCALE = str(METAPLOT_CFG.get("axis_scale", "median_length"))
METAPLOT_AXIS_WIDTHS = METAPLOT_CFG.get("axis_widths", [1, 1, 1])
# By default preserve the experimental IP/Input comparisons. An explicit groups
# mapping has the same format as bigwigCompare.groups; filtered sample subsets
# are supported without pooling unrelated conditions into an all-sample mean.
METAPLOT_GROUPS = METAPLOT_CFG.get("groups", config.get("bigwigCompare", {}).get("groups", {}))
METAPLOT_TARGETS = []
if METAPLOT_ENABLED:
    METAPLOT_TARGETS = expand(
        f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.{METAPLOT_FORMAT}",
        sample=METAPLOT_SAMPLES,
    ) + [
        f"{OUTDIR}/metaplot/integrated.UTR_CDS_UTR.{METAPLOT_FORMAT}",
        f"{OUTDIR}/metaplot/integrated.UTR_CDS_UTR.profile.tsv",
        f"{OUTDIR}/metaplot/reference.selection.tsv",
        f"{OUTDIR}/metaplot/reference.summary.json",
    ] + expand(
        f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.stats.json",
        sample=METAPLOT_SAMPLES,
    ) + expand(
        f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.gene_signal.tsv",
        sample=METAPLOT_SAMPLES,
    )


def metaplot_track_specifications():
    specifications = []
    index = 0
    for track in METAPLOT_TRACKS:
        stranded = track.startswith(("R2firstbaseCPM.", "crosslinkCPM."))
        specifications.append({"plus": index, "minus": index + int(stranded)})
        index += 2 if stranded else 1
    return specifications


def metaplot_bigwig_paths(sample):
    paths = []
    for track in METAPLOT_TRACKS:
        base = f"{OUTDIR}/bigwig/{sample}/{sample}.{track}"
        if track.startswith(("R2firstbaseCPM.", "crosslinkCPM.")):
            paths += [base + ".plus.bw", base + ".minus.bw"]
        else:
            paths.append(base + ".bw")
    return paths


def metaplot_bigwigs(wildcards):
    return metaplot_bigwig_paths(wildcards.sample)


if METAPLOT_ENABLED:
    rule prepare_metagene_regions:
        input:
            gtf=config["reference"]["gtf"],
            bws=[path for sample in METAPLOT_SAMPLES for path in metaplot_bigwig_paths(sample)],
            helper=workflow.source_path("scripts/metagene.py")
        output:
            models=f"{OUTDIR}/metaplot/reference.transcripts.jsonl",
            selection=f"{OUTDIR}/metaplot/reference.selection.tsv",
            summary=f"{OUTDIR}/metaplot/reference.summary.json"
        params:
            min_region_length=int(METAPLOT_CFG.get("min_region_length", 1))
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/gtf_to_gene_bed.py"

    rule plot_metaplot_sample:
        input:
            bws=metaplot_bigwigs,
            models=rules.prepare_metagene_regions.output.models,
            annotation_summary=rules.prepare_metagene_regions.output.summary,
            helper=workflow.source_path("scripts/metagene.py")
        output:
            plot=f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.{METAPLOT_FORMAT}",
            data=f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.profile.tsv",
            stats=f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.stats.json",
            gene_signal=f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.gene_signal.tsv"
        resources:
            mem_mb=2000
        conda:
            "envs/qc.yaml"
        params:
            tracks=METAPLOT_TRACKS,
            track_specs=metaplot_track_specifications(),
            bins=METAPLOT_BINS,
            axis_scale=METAPLOT_AXIS_SCALE,
            axis_widths=METAPLOT_AXIS_WIDTHS
        script:
            "scripts/compute_segmented_metaplot.py"

    rule plot_metaplot_integrated:
        input:
            profiles=expand(
                f"{OUTDIR}/metaplot/{{sample}}/{{sample}}.UTR_CDS_UTR.profile.tsv",
                sample=METAPLOT_SAMPLES,
            ),
            helper=workflow.source_path("scripts/metagene.py")
        output:
            plot=f"{OUTDIR}/metaplot/integrated.UTR_CDS_UTR.{METAPLOT_FORMAT}",
            data=f"{OUTDIR}/metaplot/integrated.UTR_CDS_UTR.profile.tsv"
        conda:
            "envs/qc.yaml"
        params:
            samples=METAPLOT_SAMPLES,
            tracks=METAPLOT_TRACKS,
            groups=METAPLOT_GROUPS,
            axis_scale=METAPLOT_AXIS_SCALE
        script:
            "scripts/plot_combined_metaplot.py"
