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
