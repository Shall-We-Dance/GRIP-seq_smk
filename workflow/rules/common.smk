# workflow/rules/common.smk
import os

OUTDIR = config["output"]["dir"]
SAMPLES = list(config["samples"].keys())
FASTQ_SAMPLES = [s for s in SAMPLES if not config["samples"][s].get("bam")]
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
        if not isinstance(sample, str) or not __import__("re").fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise ValueError(f"Unsafe sample name: {sample!r}")
        if sample_cfg.get("bam"):
            if not isinstance(sample_cfg["bam"], str):
                raise ValueError(f"samples.{sample}.bam must be a path string")
            if sample_cfg.get("R1") or sample_cfg.get("R2"):
                raise ValueError(f"samples.{sample}: choose BAM or FASTQ input, not both")
            continue
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

    distribution = config.get("motif_distribution", {})
    if distribution.get("enable", False):
        chosen = distribution.get("samples", "all")
        chosen = SAMPLES if chosen == "all" else chosen
        if not isinstance(chosen, list) or not chosen or set(chosen) - set(SAMPLES):
            raise ValueError("motif_distribution.samples must be all or known sample names")
        window = int(distribution.get("window", 50))
        bg_min = int(distribution.get("background_min_distance", 200))
        bg_max = int(distribution.get("background_max_distance", 1000))
        if window < 1 or bg_min <= 2 * window or bg_max < bg_min:
            raise ValueError("Motif background distances must satisfy max>=min>2*window>=2")
        seen = set()
        for motif in distribution.get("motifs", [{"name": "DRACH", "sequence": "DRACH", "anchor_pos": 3}]):
            name = motif.get("name", "")
            sequence = motif.get("sequence", "").upper()
            if not name or name in seen or not sequence or set(sequence)-set("ACGTRYSWKMBDHVN"):
                raise ValueError("Motifs need distinct names and valid DNA IUPAC sequences")
            if not 1 <= int(motif.get("anchor_pos", 1)) <= len(sequence):
                raise ValueError("Motif anchor_pos must be within its sequence")
            seen.add(name)
    if dmc.get("enable", False) or distribution.get("enable", False):
        if not 0 < float(dmc.get("holdout_fraction", 0.2)) < 1:
            raise ValueError("de_novo_motif.holdout_fraction must be between zero and one")
        if int(dmc.get("max_sequences", 2000)) < int(dmc.get("min_sequences", 10)):
            raise ValueError("de_novo_motif.max_sequences must be >= min_sequences")
        if dmc.get("sequence_center", "r2_start") not in ("r2_start", "crosslink"):
            raise ValueError("de_novo_motif.sequence_center must be r2_start or crosslink")
        if not 1 <= int(dmc.get("min_width", 4)) <= int(dmc.get("max_width", 8)) <= 2*int(dmc.get("sequence_flank", 10))+1:
            raise ValueError("MEME motif widths must fit the extracted sequence length")

    metaplot = config.get("metaplot", {})
    if metaplot.get("enable", False):
        if not config.get("reference", {}).get("gtf"):
            raise ValueError("metaplot is enabled but reference.gtf is empty.")
        metaplot_samples = metaplot.get("samples", "all")
        metaplot_samples = SAMPLES if metaplot_samples == "all" else metaplot_samples
        unknown = sorted(set(metaplot_samples) - set(SAMPLES))
        if unknown:
            raise ValueError(f"metaplot.samples contains unknown samples: {unknown}")
        allowed = {f"{track}.{variant}" for track in ("bamCPM", "R2firstbaseCPM", "crosslinkCPM") for variant in ("noblacklist", "blacklist")}
        tracks = metaplot.get("tracks", ["R2firstbaseCPM.noblacklist", "crosslinkCPM.noblacklist"])
        if not isinstance(tracks, list) or not tracks or len(set(tracks)) != len(tracks) or set(tracks) - allowed:
            raise ValueError("metaplot.tracks must be a nonempty list of supported, distinct signal tracks")
        if not FILTER_BLACKLIST and any(t.endswith(".blacklist") for t in tracks):
            raise ValueError("Blacklist metaplot tracks require filter_blacklist=true")
        if metaplot.get("axis_scale", "median_length") not in ("median_length", "fixed"):
            raise ValueError("metaplot.axis_scale must be median_length or fixed")
        for segment in ("utr5", "cds", "utr3"):
            if int(metaplot.get("bins", {}).get(segment, 100)) < 1:
                raise ValueError(f"metaplot.bins.{segment} must be positive")
        if int(metaplot.get("min_region_length", 1)) < 1:
            raise ValueError("metaplot.min_region_length must be positive")

    signal = config.get("signal", {})
    if signal.get("five_prime_clip", "exclude") not in ("exclude", "aligned"):
        raise ValueError("signal.five_prime_clip must be exclude or aligned")
    if int(signal.get("crosslink_shift", 1)) < 0:
        raise ValueError("signal.crosslink_shift must be nonnegative")
    if not 0 <= int(config.get("filtering", {}).get("min_mapq", 11)) <= 255:
        raise ValueError("filtering.min_mapq must be between 0 and 255")
    for key, value in config.get("threads", {}).items():
        if int(value) < 1:
            raise ValueError(f"threads.{key} must be positive")
    if bwc.get("enable", True):
        if int(bwc.get("binSize", 1)) < 1 or float(bwc.get("pseudocount", 1)) <= 0:
            raise ValueError("bigwigCompare.binSize and pseudocount must be positive")
        if bwc.get("operation", "log2") not in ("log2", "ratio", "subtract", "add", "mean", "reciprocal_ratio", "first", "second"):
            raise ValueError("Unsupported bigwigCompare.operation")
        if not FILTER_BLACKLIST and any(t.endswith(".blacklist") for t in bwc.get("tracks", [])):
            raise ValueError("Blacklist comparison tracks require filter_blacklist=true")
        for name in bwc.get("groups", {}):
            if not __import__("re").fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
                raise ValueError(f"Unsafe comparison group name: {name!r}")


rule faidx_reference:
    input:
        fa=config["reference"]["fasta"]
    output:
        fai=config["reference"]["fasta"] + ".fai"
    threads: 2
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input.fa:q}"

validate_workflow_config()

if FILTER_BLACKLIST and blacklist_config().get("url"):
    rule get_blacklist:
        output:
            bed=blacklist_path()
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/get_blacklist.py"
