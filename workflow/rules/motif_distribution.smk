# motif_de_novo.smk declares the sequence-unbiased RT calling rule and shared
# discovery/holdout settings. Include it before this file.
MOTIF_DISTRIBUTION_TARGETS = []
if MOTIF_DIST_ENABLED:
    MOTIF_DISTRIBUTION_TARGETS = [
        f"{OUTDIR}/motif_distribution/{sample}/{sample}.{suffix}"
        for sample in MOTIF_DIST_SAMPLES
        for suffix in ("offsets.tsv", "r2_offsets.tsv", "statistics.json", "distribution.png", "distribution.pdf", "background_pairs.tsv")
    ]

    rule motif_offset_distribution:
        input:
            bed=f"{OUTDIR}/motif_de_novo/{{sample}}/{{sample}}.rt_termination.bed",
            fasta=config["reference"]["fasta"],
            helper=os.path.join(MOTIF_SCRIPT_DIR, "motif_utils.py"),
            fai=rules.faidx_reference.output.fai,
            blacklist=[blacklist_path()] if MOTIF_VARIANT == "blacklist" else [],
            meme=lambda wc: [f"{OUTDIR}/motif_de_novo/{wc.sample}/{wc.sample}.meme.txt"]
                if DE_NOVO_ENABLED and bool(MOTIF_DIST_CFG.get("scan_de_novo", True)) and wc.sample in MOTIF_SAMPLES else []
        output:
            table=f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.offsets.tsv",
            r2_table=f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.r2_offsets.tsv",
            stats=f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.statistics.json",
            png=f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.distribution.png",
            pdf=f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.distribution.pdf",
            pairs=f"{OUTDIR}/motif_distribution/{{sample}}/{{sample}}.background_pairs.tsv"
        log:
            f"{OUTDIR}/logs/motif_distribution/{{sample}}.log"
        conda:
            "envs/qc.yaml"
        threads: 1
        resources:
            mem_mb=4000
        params:
            script_dir=MOTIF_SCRIPT_DIR,
            variant=MOTIF_VARIANT,
            motifs=MOTIF_DIST_CFG.get("motifs", [{"name": "DRACH", "sequence": "DRACH", "anchor_pos": 3}]),
            window=int(MOTIF_DIST_CFG.get("window", 50)),
            site_shift=int(config.get("signal", {}).get("crosslink_shift", 1)),
            seed=int(DE_NOVO_CFG.get("seed", 17)),
            holdout_fraction=float(DE_NOVO_CFG.get("holdout_fraction", 0.2)),
            holdout_block_size=int(DE_NOVO_CFG.get("holdout_block_size", 10000)),
            holdout_guard=int(DE_NOVO_CFG.get("holdout_guard", 100)),
            consensus_mass=float(MOTIF_DIST_CFG.get("consensus_probability_mass", 0.8)),
            background_min_distance=int(MOTIF_DIST_CFG.get("background_min_distance", 200)),
            background_max_distance=int(MOTIF_DIST_CFG.get("background_max_distance", 1000)),
            background_attempts=int(MOTIF_DIST_CFG.get("background_attempts", 20))
        script:
            "scripts/motif_distribution.py"
