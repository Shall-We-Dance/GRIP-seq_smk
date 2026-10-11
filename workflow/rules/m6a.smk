"""Optional Input-enriched, motif-supported m6A candidates and set comparisons."""
import math
import re

M6A_CFG = config.get("m6a", {})
M6A_ENABLED = bool(M6A_CFG.get("enable", False))
M6A_TARGETS = []
if M6A_ENABLED:
    M6A_CONTROLS = M6A_CFG.get("controls", {})
    if not isinstance(M6A_CONTROLS, dict) or not M6A_CONTROLS:
        raise ValueError("m6a.controls must map each IP sample to one Input sample")
    if any(sample not in SAMPLES or control not in SAMPLES or sample == control for sample, control in M6A_CONTROLS.items()):
        raise ValueError("m6a.controls requires distinct known IP and Input samples")
    if set(M6A_CONTROLS) & set(M6A_CONTROLS.values()):
        raise ValueError("An m6a sample cannot have both IP and Input roles")
    M6A_SAMPLES = list(M6A_CONTROLS)
    M6A_GROUPS = M6A_CFG.get("groups", {})
    if not isinstance(M6A_GROUPS, dict) or not M6A_GROUPS:
        raise ValueError("m6a.groups must define at least one group of IP libraries")
    for name, members in M6A_GROUPS.items():
        if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", str(name)) or name == "Atlas" or
                not isinstance(members, list) or not members or len(set(members)) != len(members) or
                set(members)-set(M6A_SAMPLES)):
            raise ValueError("m6a.groups needs safe names (other than Atlas) and distinct known IP libraries")
    if set(M6A_CFG.get("exclude_samples", []))-set(M6A_SAMPLES):
        raise ValueError("m6a.exclude_samples contains an unknown IP library")
    for key, default, lower in (("min_reads", 10, 0), ("min_library_support", 2, 1),
                                ("metagene_bins", 30, 1), ("offset_window", 25, 1)):
        if int(M6A_CFG.get(key, default)) < lower:
            raise ValueError(f"m6a.{key} must be >= {lower}")
    for key, default, lower in (("neighbor_fold", 1.5, 1), ("min_input_fold", 2, 0), ("pseudocount", 0.5, 0)):
        value = float(M6A_CFG.get(key, default))
        if not math.isfinite(value) or value <= lower:
            raise ValueError(f"m6a.{key} must be finite and > {lower}")
    int(M6A_CFG.get("adenosine_offset", 3))
    M6A_VARIANT = M6A_CFG.get("variant", "auto")
    if M6A_VARIANT == "auto":
        M6A_VARIANT = "blacklist" if FILTER_BLACKLIST else "noblacklist"
    if M6A_VARIANT not in SIGNAL_VARIANTS:
        raise ValueError("m6a.variant is not an enabled signal variant")
    M6A_SCRIPT_DIR = os.path.join(workflow.basedir, "rules", "scripts")
    M6A_TARGETS = [f"{OUTDIR}/m6a/comparison"]

    rule call_m6a_candidates:
        input:
            counts=f"{OUTDIR}/signal/{{sample}}/{{sample}}.R2firstbase.{M6A_VARIANT}.counts.bed.gz",
            mapping=f"{OUTDIR}/signal/{{sample}}/{{sample}}.R2_to_crosslink.{M6A_VARIANT}.counts.tsv.gz",
            crosslinks=f"{OUTDIR}/signal/{{sample}}/{{sample}}.crosslink.{M6A_VARIANT}.counts.bed.gz",
            control=lambda wc: f"{OUTDIR}/signal/{M6A_CONTROLS[wc.sample]}/{M6A_CONTROLS[wc.sample]}.crosslink.{M6A_VARIANT}.counts.bed.gz",
            ip_stats=f"{OUTDIR}/signal/{{sample}}/{{sample}}.noblacklist.signal.stats.json",
            control_stats=lambda wc: f"{OUTDIR}/signal/{M6A_CONTROLS[wc.sample]}/{M6A_CONTROLS[wc.sample]}.noblacklist.signal.stats.json",
            fasta=config["reference"]["fasta"],
            fai=rules.faidx_reference.output.fai,
            helpers=[os.path.join(M6A_SCRIPT_DIR, name) for name in ("call_rt_termination_sites.py", "motif_utils.py")]
        output:
            candidates=f"{OUTDIR}/m6a/{{sample}}/candidates.tsv",
            bed=f"{OUTDIR}/m6a/{{sample}}/candidates.bed",
            all_rt=f"{OUTDIR}/m6a/{{sample}}/all_RT_sites.tsv",
            offsets=f"{OUTDIR}/m6a/{{sample}}/DRACH_offsets.tsv",
            stats=f"{OUTDIR}/m6a/{{sample}}/stats.json"
        wildcard_constraints:
            sample="|".join(re.escape(sample) for sample in M6A_SAMPLES)
        params:
            settings=M6A_CFG,
            script_dir=M6A_SCRIPT_DIR,
            variant=M6A_VARIANT,
            control=lambda wc: M6A_CONTROLS[wc.sample],
            crosslink_shift=int(config.get("signal", {}).get("crosslink_shift", 1))
        conda:
            "envs/py_signal.yaml"
        resources:
            mem_mb=4000
        script:
            "scripts/m6a_sites.py"

    rule compare_m6a_candidates:
        input:
            calls=expand(f"{OUTDIR}/m6a/{{sample}}/candidates.tsv", sample=M6A_SAMPLES),
            offsets=expand(f"{OUTDIR}/m6a/{{sample}}/DRACH_offsets.tsv", sample=M6A_SAMPLES),
            stats=expand(f"{OUTDIR}/m6a/{{sample}}/stats.json", sample=M6A_SAMPLES),
            fasta=config["reference"]["fasta"],
            fai=rules.faidx_reference.output.fai,
            gtf=[config["reference"]["gtf"]] if config["reference"].get("gtf") else [],
            atlas=[M6A_CFG["atlas_bed"]] if M6A_CFG.get("atlas_bed") else [],
            helpers=[os.path.join(M6A_SCRIPT_DIR, name) for name in ("m6a_sites.py", "metagene.py", "call_rt_termination_sites.py", "motif_utils.py")]
        output:
            report=directory(f"{OUTDIR}/m6a/comparison")
        params:
            settings=M6A_CFG,
            samples=M6A_SAMPLES,
            script_dir=M6A_SCRIPT_DIR
        conda:
            "envs/m6a.yaml"
        resources:
            mem_mb=6000
        script:
            "scripts/compare_m6a_sites.py"
