# workflow/rules/tracks.smk
import os

OUTDIR = config["output"]["dir"]

if FILTER_BLACKLIST:
    rule make_tracks:
        input:
            bw1=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.noblacklist.bw",
            bw2=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.blacklist.bw",
            bw3=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.noblacklist.bw",
            bw4=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.blacklist.bw",
        output:
            txt=f"{OUTDIR}/tracks/{{sample}}/{{sample}}.tracks.txt"
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/make_tracks.py"
else:
    rule make_tracks:
        input:
            bw1=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.bamCPM.noblacklist.bw",
            bw3=f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.R2firstbaseCPM.noblacklist.bw",
        output:
            txt=f"{OUTDIR}/tracks/{{sample}}/{{sample}}.tracks.txt"
        conda:
            "envs/py_signal.yaml"
        script:
            "scripts/make_tracks.py"
