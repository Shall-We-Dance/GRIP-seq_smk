# workflow/rules/common.smk
import os

OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())
FILTER_BLACKLIST = bool(config.get("filter_blacklist", False))


def blacklist_config():
    return config.get("blacklist", {})


def blacklist_path():
    blacklist = blacklist_config()
    if blacklist.get("path"):
        return blacklist["path"]
    if blacklist.get("url"):
        cache_dir = blacklist.get("cache_dir", "resources/blacklist")
        filename = os.path.basename(blacklist["url"])
        if filename.endswith(".gz"):
            filename = filename[:-3]
        return os.path.join(cache_dir, filename)
    return None


def validate_workflow_config():
    if not SAMPLES:
        raise ValueError("config.samples must contain at least one sample.")

    for sample, sample_cfg in config.get("samples", {}).items():
        r1 = sample_cfg.get("R1", [])
        r2 = sample_cfg.get("R2", [])
        if not isinstance(r1, list) or not r1:
            raise ValueError(f"samples.{sample}.R1 must be a non-empty list.")
        if not isinstance(r2, list) or not r2:
            raise ValueError(f"samples.{sample}.R2 must be a non-empty list.")
        if len(r1) != len(r2):
            raise ValueError(
                f"samples.{sample}.R1 and R2 must have the same number of FASTQs."
            )

    blacklist = blacklist_config()
    has_blacklist_path = bool(blacklist.get("path"))
    has_blacklist_url = bool(blacklist.get("url"))
    if FILTER_BLACKLIST and has_blacklist_path == has_blacklist_url:
        raise ValueError(
            "filter_blacklist is enabled; provide exactly one of "
            "blacklist.path or blacklist.url."
        )

    bwc = config.get("bigwigCompare", {})
    if bwc.get("enable", True):
        groups = bwc.get("groups", {})
        if not groups:
            raise ValueError(
                "bigwigCompare is enabled; define at least one comparison in "
                "bigwigCompare.groups."
            )
        for gname, g in groups.items():
            if not isinstance(g, dict):
                raise ValueError(
                    f"bigwigCompare.groups.{gname} must be a dict with IP/Input lists."
                )
            for key in ("IP", "Input"):
                lst = g.get(key, [])
                if not isinstance(lst, list) or not lst:
                    raise ValueError(
                        f"bigwigCompare.groups.{gname}.{key} must be a non-empty "
                        "list of sample names."
                    )
                for s in lst:
                    if s not in config["samples"]:
                        raise ValueError(
                            f"bigwigCompare.groups.{gname}.{key}: sample '{s}' "
                            "not found in config.samples."
                        )

    mc = config.get("motif_anchoring", {})
    if mc.get("enable", True):
        motif = str(mc.get("motif", "DRACH"))
        motif_pos = int(mc.get("motif_pos", 3))
        if not 1 <= motif_pos <= len(motif):
            raise ValueError(
                f"motif_anchoring.motif_pos ({motif_pos}) must be within "
                f"1..{len(motif)} of motif '{motif}'."
            )
        if int(mc.get("search_window", 0)) < 0:
            raise ValueError("motif_anchoring.search_window must be >= 0.")

    dmc = config.get("de_novo_motif", {})
    if dmc.get("enable", False):
        motif_samples = dmc.get("samples", "all")
        motif_samples = SAMPLES if motif_samples == "all" else motif_samples
        unknown = sorted(set(motif_samples) - set(SAMPLES))
        if unknown:
            raise ValueError(f"de_novo_motif.samples contains unknown samples: {unknown}")
        if int(dmc.get("min_reads", 40)) < 1:
            raise ValueError("de_novo_motif.min_reads must be >= 1.")
        if float(dmc.get("neighbor_fold", 1.5)) <= 1:
            raise ValueError("de_novo_motif.neighbor_fold must be > 1.")

    metaplot = config.get("metaplot", {})
    if metaplot.get("enable", False):
        if not config.get("reference", {}).get("gtf"):
            raise ValueError("metaplot is enabled but reference.gtf is empty.")
        metaplot_samples = metaplot.get("samples", "all")
        metaplot_samples = SAMPLES if metaplot_samples == "all" else metaplot_samples
        unknown = sorted(set(metaplot_samples) - set(SAMPLES))
        if unknown:
            raise ValueError(f"metaplot.samples contains unknown samples: {unknown}")
        tracks = metaplot.get("tracks", [])
        required_tracks = {"bamCPM.noblacklist", "R2firstbaseCPM.noblacklist"}
        if set(tracks) != required_tracks:
            raise ValueError(
                "metaplot.tracks must contain bamCPM.noblacklist and "
                "R2firstbaseCPM.noblacklist."
            )
        for key in ("upstream", "downstream", "gene_body_length", "bin_size"):
            if int(metaplot.get(key, 1)) < 1:
                raise ValueError(f"metaplot.{key} must be >= 1.")

rule faidx_reference:
    input:
        fa=config["reference"]["fasta"]
    output:
        fai=config["reference"]["fasta"] + ".fai"
    threads: 2
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input.fa}"

validate_workflow_config()

if FILTER_BLACKLIST and blacklist_config().get("url"):
    rule get_blacklist:
        output:
            bed=blacklist_path()
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/get_blacklist.py"
