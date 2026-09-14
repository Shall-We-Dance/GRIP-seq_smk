import os
import re

OUTDIR = config["output"]["dir"]
DE_NOVO_CFG = config.get("de_novo_motif", {})
DE_NOVO_ENABLED = bool(DE_NOVO_CFG.get("enable", False))
MOTIF_DIST_CFG = config.get("motif_distribution", {})
MOTIF_DIST_ENABLED = bool(MOTIF_DIST_CFG.get("enable", False))
MOTIF_VARIANT = str(DE_NOVO_CFG.get("variant", "auto"))
if MOTIF_VARIANT == "auto":
    MOTIF_VARIANT = "blacklist" if FILTER_BLACKLIST else "noblacklist"
if (DE_NOVO_ENABLED or MOTIF_DIST_ENABLED) and (MOTIF_VARIANT not in ("noblacklist", "blacklist") or (MOTIF_VARIANT == "blacklist" and not FILTER_BLACKLIST)):
    raise ValueError("de_novo_motif.variant must exist in the enabled signal variants")
MOTIF_SCRIPT_DIR = os.path.join(workflow.basedir, "rules", "scripts")
configured_motif_samples = DE_NOVO_CFG.get("samples", "all")
MOTIF_SAMPLES = SAMPLES if configured_motif_samples == "all" else list(configured_motif_samples)
configured_distribution_samples = MOTIF_DIST_CFG.get("samples", "all")
MOTIF_DIST_SAMPLES = SAMPLES if configured_distribution_samples == "all" else list(configured_distribution_samples)

# Derive pools from explicit biological roles, never from substrings in names.
# Controls shared across comparisons are included only once.
comparison_groups = config.get("bigwigCompare", {}).get("groups", {})
default_pools = {}
for side in ("IP", "Input"):
    members = sorted({sample for group in comparison_groups.values() for sample in group.get(side, []) if sample in MOTIF_SAMPLES})
    if members:
        default_pools[side] = members
MOTIF_POOLS = DE_NOVO_CFG.get("pools", default_pools)
for pool, members in MOTIF_POOLS.items():
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(pool)) or not isinstance(members, list) or not members:
        raise ValueError("de_novo_motif.pools needs safe pool names and nonempty sample lists")
    if set(members) - set(MOTIF_SAMPLES):
        raise ValueError(f"Unknown/disabled motif samples in pool {pool}")
known_ip = {s for group in comparison_groups.values() for s in group.get("IP", [])}
known_input = {s for group in comparison_groups.values() for s in group.get("Input", [])}
if known_ip & known_input:
    raise ValueError("A sample cannot be both IP and Input for motif analysis")
for pool, members in MOTIF_POOLS.items():
    if set(members) & known_ip and set(members) & known_input:
        raise ValueError(f"Motif pool {pool} mixes IP and Input samples")

DE_NOVO_MOTIF_TARGETS = []
if DE_NOVO_ENABLED:
    DE_NOVO_MOTIF_TARGETS = [
        f"{OUTDIR}/motif_de_novo/{sample}/{sample}.{suffix}"
        for sample in MOTIF_SAMPLES
        for suffix in ("meme.html", "meme.txt", "motif1.png", "meme.stats.json")
    ] + [
        f"{OUTDIR}/motif_de_novo/pooled/{pool}/{pool}.{suffix}"
        for pool in MOTIF_POOLS
        for suffix in ("meme.html", "meme.txt", "motif1.png", "meme.stats.json")
    ]

if DE_NOVO_ENABLED or MOTIF_DIST_ENABLED:
    rule call_rt_termination_sites:
        input:
            counts=f"{OUTDIR}/signal/{{sample}}/{{sample}}.R2firstbase.{MOTIF_VARIANT}.counts.bed.gz",
            mapping=f"{OUTDIR}/signal/{{sample}}/{{sample}}.R2_to_crosslink.{MOTIF_VARIANT}.counts.tsv.gz",
            fasta=config["reference"]["fasta"],
            helper=os.path.join(MOTIF_SCRIPT_DIR, "motif_utils.py"),
            fai=rules.faidx_reference.output.fai
        output:
            bed=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.bed",
            fasta=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.fa",
            stats=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.stats.tsv"
        log:
            f"{OUTDIR}/logs/motif_de_novo/{{sample}}.sites.log"
        conda:
            "envs/py_signal.yaml"
        threads: 1
        resources:
            mem_mb=8000
        params:
            script_dir=MOTIF_SCRIPT_DIR,
            variant=MOTIF_VARIANT,
            min_reads=int(DE_NOVO_CFG.get("min_reads", 40)),
            neighbor_fold=float(DE_NOVO_CFG.get("neighbor_fold", 1.5)),
            flank=int(DE_NOVO_CFG.get("sequence_flank", 10)),
            sequence_center=str(DE_NOVO_CFG.get("sequence_center", "r2_start")),
            site_shift=int(config.get("signal", {}).get("crosslink_shift", 1))
        script:
            "scripts/call_rt_termination_sites.py"

if DE_NOVO_ENABLED:
    rule meme_de_novo_sample:
        input:
            fasta=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.fa",
            helper=os.path.join(MOTIF_SCRIPT_DIR, "motif_utils.py")
        output:
            html=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.meme.html",
            txt=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.meme.txt",
            logo=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.motif1.png",
            sequences=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.discovery.fa",
            stats=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.meme.stats.json"
        log:
            f"{OUTDIR}/logs/motif_de_novo/{{sample}}.meme.log"
        conda:
            "envs/motif.yaml"
        threads: int(DE_NOVO_CFG.get("threads", 4))
        resources:
            mem_mb=4000
        params:
            script_dir=MOTIF_SCRIPT_DIR,
            label=lambda wc: wc.sample,
            min_sequences=int(DE_NOVO_CFG.get("min_sequences", 10)),
            max_sequences=int(DE_NOVO_CFG.get("max_sequences", 2000)),
            holdout_fraction=float(DE_NOVO_CFG.get("holdout_fraction", 0.2)),
            holdout_block_size=int(DE_NOVO_CFG.get("holdout_block_size", 10000)),
            holdout_guard=int(DE_NOVO_CFG.get("holdout_guard", 100)),
            seed=int(DE_NOVO_CFG.get("seed", 17)),
            nmotifs=int(DE_NOVO_CFG.get("nmotifs", 5)),
            min_width=int(DE_NOVO_CFG.get("min_width", 4)),
            max_width=int(DE_NOVO_CFG.get("max_width", 8))
        script:
            "scripts/run_meme.py"

    rule meme_de_novo_pool:
        input:
            fasta=lambda wc: [f"{OUTDIR}/motif_de_novo/{sample}/{sample}.rt_termination.fa" for sample in MOTIF_POOLS[wc.pool]],
            helper=os.path.join(MOTIF_SCRIPT_DIR, "motif_utils.py")
        output:
            html=f"{OUTDIR}/motif_de_novo/pooled/{{pool}}/{{pool}}.meme.html",
            txt=f"{OUTDIR}/motif_de_novo/pooled/{{pool}}/{{pool}}.meme.txt",
            logo=f"{OUTDIR}/motif_de_novo/pooled/{{pool}}/{{pool}}.motif1.png",
            sequences=f"{OUTDIR}/motif_de_novo/pooled/{{pool}}/{{pool}}.discovery.fa",
            stats=f"{OUTDIR}/motif_de_novo/pooled/{{pool}}/{{pool}}.meme.stats.json"
        wildcard_constraints:
            pool="|".join(re.escape(pool) for pool in MOTIF_POOLS) or "(?!)"
        log:
            f"{OUTDIR}/logs/motif_de_novo/pooled/{{pool}}.meme.log"
        conda:
            "envs/motif.yaml"
        threads: int(DE_NOVO_CFG.get("threads", 4))
        resources:
            mem_mb=4000
        params:
            script_dir=MOTIF_SCRIPT_DIR,
            label=lambda wc: f"Pooled {wc.pool} (unique genomic loci)",
            min_sequences=int(DE_NOVO_CFG.get("min_sequences", 10)),
            max_sequences=int(DE_NOVO_CFG.get("max_sequences", 2000)),
            holdout_fraction=float(DE_NOVO_CFG.get("holdout_fraction", 0.2)),
            holdout_block_size=int(DE_NOVO_CFG.get("holdout_block_size", 10000)),
            holdout_guard=int(DE_NOVO_CFG.get("holdout_guard", 100)),
            seed=int(DE_NOVO_CFG.get("seed", 17)),
            nmotifs=int(DE_NOVO_CFG.get("nmotifs", 5)),
            min_width=int(DE_NOVO_CFG.get("min_width", 4)),
            max_width=int(DE_NOVO_CFG.get("max_width", 8))
        script:
            "scripts/run_meme.py"
