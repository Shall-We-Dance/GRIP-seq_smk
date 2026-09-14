# BAM coverage can reuse externally validated baselines during BAM reanalysis.
import os

OUTDIR = config["output"]["dir"]

for _sample, _settings in config.get("samples", {}).items():
    if "bam_coverage" not in _settings:
        continue
    if not _settings.get("bam"):
        raise ValueError(f"samples.{_sample}.bam_coverage requires BAM input mode")
    _coverage = _settings["bam_coverage"]
    if not isinstance(_coverage, dict) or set(_coverage) - {"noblacklist", "blacklist"}:
        raise ValueError(f"samples.{_sample}.bam_coverage must map noblacklist/blacklist to source BigWigs")
    if any(not isinstance(path, str) or not path.strip() for path in _coverage.values()):
        raise ValueError(f"samples.{_sample}.bam_coverage paths must be nonempty strings")


def coverage_reuse_input(sample, variant):
    path = config.get("samples", {}).get(sample, {}).get("bam_coverage", {}).get(variant)
    if path:
        destination = f"{OUTDIR}/bigwig/{sample}/{sample}.bamCPM.{variant}.bw"
        if os.path.realpath(path) == os.path.realpath(destination):
            raise ValueError(f"samples.{sample}.bam_coverage.{variant} must differ from its output")
    return [path] if path else []


rule bam_to_bigwig_cpm_noblacklist:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
        reuse=lambda wc: coverage_reuse_input(wc.sample, "noblacklist")
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.noblacklist.bw"
    log:
        f"{OUTDIR}/logs/bigwig/{{sample}}.bamCoverage.noblacklist.log"
    threads: lambda wc: 1 if coverage_reuse_input(wc.sample, "noblacklist") else int(config.get("threads", {}).get("bigwig", 8))
    conda:
        "envs/qc.yaml"
    script:
        "scripts/bam_coverage.py"

if FILTER_BLACKLIST:
    rule bam_to_bigwig_cpm_blacklist:
        input:
            bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
            bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
            bl=blacklist_path(),
            reuse=lambda wc: coverage_reuse_input(wc.sample, "blacklist")
        output:
            bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.blacklist.bw"
        log:
            f"{OUTDIR}/logs/bigwig/{{sample}}.bamCoverage.blacklist.log"
        threads: lambda wc: 1 if coverage_reuse_input(wc.sample, "blacklist") else int(config.get("threads", {}).get("bigwig", 8))
        conda:
            "envs/qc.yaml"
        script:
            "scripts/bam_coverage.py"

# A single BAM traversal creates endpoints, inferred cross-links and audit counts.
SIGNAL_CFG = config.get("signal", {})
SIGNAL_VARIANTS = ["noblacklist"] + (["blacklist"] if FILTER_BLACKLIST else [])
SIGNAL_TARGETS = expand(
    f"{OUTDIR}/signal/{{sample}}/{{sample}}.{{variant}}.signal.stats.json",
    sample=SAMPLES, variant=SIGNAL_VARIANTS,
)

rule r2_first_base_signal:
    input:
        bam=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam",
        bai=f"{OUTDIR}/star/{{sample}}/{{sample}}.unique.mapq11.sorted.bam.bai",
        fai=rules.faidx_reference.output.fai,
        bl=lambda wc: [blacklist_path()] if wc.variant == "blacklist" else []
    output:
        bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.{{variant}}.bw",
        plus=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.{{variant}}.plus.bw",
        minus=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.{{variant}}.minus.bw",
        crosslink_bw=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.crosslinkCPM.{{variant}}.bw",
        crosslink_plus=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.crosslinkCPM.{{variant}}.plus.bw",
        crosslink_minus=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.crosslinkCPM.{{variant}}.minus.bw",
        bed=f"{OUTDIR}/signal/{{sample}}/{{sample}}.R2firstbase.{{variant}}.counts.bed.gz",
        crosslink_bed=f"{OUTDIR}/signal/{{sample}}/{{sample}}.crosslink.{{variant}}.counts.bed.gz",
        mapping=f"{OUTDIR}/signal/{{sample}}/{{sample}}.R2_to_crosslink.{{variant}}.counts.tsv.gz",
        stats=f"{OUTDIR}/signal/{{sample}}/{{sample}}.{{variant}}.signal.stats.json"
    wildcard_constraints:
        variant="|".join(SIGNAL_VARIANTS)
    log:
        f"{OUTDIR}/logs/5prime/{{sample}}.{{variant}}.signal.log"
    conda:
        "envs/py_signal.yaml"
    threads: 1
    resources:
        mem_mb=2048
    params:
        min_mapq=int(config["filtering"]["min_mapq"]),
        five_prime_clip=str(SIGNAL_CFG.get("five_prime_clip", "exclude")),
        crosslink_shift=int(SIGNAL_CFG.get("crosslink_shift", 1)),
        require_proper_pair=bool(SIGNAL_CFG.get("require_proper_pair", False)),
        exclude_duplicates=bool(SIGNAL_CFG.get("exclude_duplicates", False)),
        exclude_qcfail=bool(SIGNAL_CFG.get("exclude_qcfail", True))
    script:
        "scripts/extract_r2_first_base.py"
