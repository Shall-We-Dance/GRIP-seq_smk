OUTDIR = config["output"]["dir"]

def browser_tracks(wc):
    paths = []
    for variant in SIGNAL_VARIANTS:
        paths.append(f"{OUTDIR}/bigwig/{wc.sample}/{wc.sample}.bamCPM.{variant}.bw")
        for signal in ("R2firstbaseCPM", "crosslinkCPM"):
            for strand in ("", ".plus", ".minus"):
                paths.append(f"{OUTDIR}/bigwig/{wc.sample}/{wc.sample}.{signal}.{variant}{strand}.bw")
    return paths

rule make_tracks:
    input:
        tracks=browser_tracks
    output:
        txt=f"{OUTDIR}/tracks/{{sample}}/{{sample}}.tracks.txt"
    conda:
        "envs/py_signal.yaml"
    script:
        "scripts/make_tracks.py"
