# workflow/rules/bigwig_compare.smk
import os
import re

OUTDIR = config["output"]["dir"]
BWC_CFG = config.get("bigwigCompare", {})
BWC_ENABLED = bool(BWC_CFG.get("enable", True))
BWC_OPERATION = str(BWC_CFG.get("operation", "log2"))
BWC_TRACKS = list(BWC_CFG.get("tracks", ["R2firstbaseCPM.noblacklist"]))
BWC_GROUPS = BWC_CFG.get("groups", {})


def bigwig_compare_targets():
    return [
        f"{OUTDIR}/bigwig/compare/{group}.{track}.{BWC_OPERATION}.bw"
        for group in BWC_GROUPS
        for track in BWC_TRACKS
    ]


BIGWIG_COMPARE_TARGETS = bigwig_compare_targets()


if BWC_ENABLED:
    # Pool (average) per-sample BigWigs of one group side (IP or Input)
    rule pool_bigwigs:
        wildcard_constraints:
            group="|".join(re.escape(g) for g in BWC_GROUPS),
            track="|".join(re.escape(t) for t in BWC_TRACKS),
            pool="IP|Input",
        input:
            bws=lambda wc: [
                f"{OUTDIR}/bigwig/{s}/{s}.{wc.track}.bw"
                for s in BWC_GROUPS[wc.group][wc.pool]
            ]
        output:
            bw=temp(f"{OUTDIR}/tmp/bigwig_compare/{{group}}.{{track}}.{{pool}}.pooled.bw")
        log:
            f"logs/bigwig_compare/{{group}}.{{track}}.{{pool}}.pool.log"
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/average_bigwigs.py"

    # log2(IP / Input) track per group
    rule bigwig_compare_log2:
        wildcard_constraints:
            group="|".join(re.escape(g) for g in BWC_GROUPS),
            track="|".join(re.escape(t) for t in BWC_TRACKS),
        input:
            bw1=f"{OUTDIR}/tmp/bigwig_compare/{{group}}.{{track}}.IP.pooled.bw",
            bw2=f"{OUTDIR}/tmp/bigwig_compare/{{group}}.{{track}}.Input.pooled.bw"
        output:
            bw=f"{OUTDIR}/bigwig/compare/{{group}}.{{track}}.{BWC_OPERATION}.bw"
        log:
            f"logs/bigwig_compare/{{group}}.{{track}}.{BWC_OPERATION}.log"
        threads: int(config.get("threads", {}).get("bigwig", 8))
        conda:
            "envs/qc.yaml"
        params:
            operation=BWC_OPERATION,
            pseudocount=float(BWC_CFG.get("pseudocount", 1.0)),
            binsize=int(BWC_CFG.get("binSize", 1))
        shell:
            r"""
            mkdir -p $(dirname {output.bw}) $(dirname {log})
            bigwigCompare \
              --bigwig1 {input.bw1} \
              --bigwig2 {input.bw2} \
              --operation {params.operation} \
              --pseudocount {params.pseudocount} \
              --binSize {params.binsize} \
              --numberOfProcessors {threads} \
              -o {output.bw} \
              > {log} 2>&1
            """
