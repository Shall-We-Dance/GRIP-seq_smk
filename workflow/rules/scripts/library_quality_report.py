"""Report endpoint precision, library depth, and expression-level reproducibility."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


def main(sm):
    outdir=Path(sm.config['output']['dir'])
    upstream=Path(sm.config.get('provenance',{}).get('upstream_results',str(outdir)))
    limits=sm.config.get('quality',{})
    minimum=int(limits.get('min_usable_r2',100000))
    clip_limit=float(limits.get('max_five_prime_clip_fraction',0.2))
    rows=[]
    for sample,path in zip(sm.config['samples'],sm.input.signals):
        stats=json.loads(Path(path).read_text())
        total=stats.get('read2_pre_endpoint_filter',stats.get('read2_total',0))
        accepted=stats.get('accepted_read2',0)
        soft=stats.get('five_prime_softclip',0)
        clipped=stats.get('excluded_five_prime_unaligned',0)
        # In aligned compatibility mode, clipping is retained and still must be reported.
        observed_clip=stats.get('five_prime_unaligned_any',clipped)
        fraction=observed_clip/total if total else 0
        flags=[]
        if accepted<minimum: flags.append('LOW_USABLE_R2')
        if fraction>clip_limit: flags.append('HIGH_5PRIME_UNALIGNED')
        raw_pairs=None
        fastp=upstream/'qc/fastp'/sample/'merged_step1.json'
        if fastp.exists():
            raw_pairs=json.loads(fastp.read_text())['summary']['before_filtering']['total_reads']//2
        rows.append(dict(sample=sample,raw_pairs=raw_pairs,primary_mapped_r2=total,
                         accepted_r2=accepted,usable_fraction=accepted/total if total else 0,
                         five_prime_softclip=soft,five_prime_unaligned_fraction=fraction,
                         excluded_five_prime_unaligned=clipped,flags=';'.join(flags) or 'PASS'))
    target=Path(sm.output.table); target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    by_track={}
    for path in sm.input.gene_signals:
        with open(path) as f:
            for r in csv.DictReader(f,delimiter='\t'):
                by_track.setdefault(r['track'],{}).setdefault(r['sample'],{})[r['gene_id']]=float(r['signal_mass'])
    correlations=[]
    heatmaps=[]
    for track, samples in sorted(by_track.items()):
        names=list(samples)
        genes=sorted(set.intersection(*(set(s) for s in samples.values()))) if samples else []
        matrix=np.array([[samples[name][g] for g in genes] for name in names])
        # Pairwise union-of-detected genes avoids thousands of shared zero genes inflating r.
        corr=np.full((len(names),len(names)),np.nan)
        for i,left in enumerate(names):
            for j,right in enumerate(names):
                keep=(matrix[i]>0)|(matrix[j]>0)
                n=int(keep.sum())
                x,y=matrix[i,keep],matrix[j,keep]
                value=float(spearmanr(x,y).statistic) if n>=3 and np.ptp(x)>0 and np.ptp(y)>0 else float('nan')
                corr[i,j]=value
                correlations.append(dict(track=track,sample_a=left,sample_b=right,n_detected_union=n,spearman_r=value))
        heatmaps.append((track,names,corr))
    with open(sm.output.correlations,'w') as f:
        writer=csv.DictWriter(f,fieldnames=['track','sample_a','sample_b','n_detected_union','spearman_r'],delimiter='\t');writer.writeheader();writer.writerows(correlations)
    summary={'samples':rows,'thresholds':{'min_usable_r2':minimum,'max_five_prime_clip_fraction':clip_limit},
             'interpretation':'QC flags are advisory; all configured samples remain included. Counts are reads, not UMI-deduplicated molecules.',
             'correlation_definition':'Spearman over per-pair union of genes with positive signal; shared zeros excluded; descriptive only.'}
    Path(sm.output.json).write_text(json.dumps(summary,indent=2)+'\n')
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(sm.output.plot) as pdf:
        fig,axes=plt.subplots(2,1,figsize=(11,8),sharex=True)
        names=[r['sample'] for r in rows];x=np.arange(len(names))
        axes[0].bar(x,[r['primary_mapped_r2'] for r in rows],color='lightgray',label='Mapped R2 before endpoint QC')
        axes[0].bar(x,[r['accepted_r2'] for r in rows],color='#218C74',label='Usable precise endpoints')
        axes[0].set_yscale('symlog',linthresh=1);axes[0].set_ylabel('R2 alignments');axes[0].legend(fontsize=8)
        axes[1].bar(x,[r['five_prime_unaligned_fraction'] for r in rows],color='#CC614A')
        axes[1].axhline(clip_limit,color='gray',ls='--');axes[1].set_ylabel('5′ unaligned fraction')
        axes[1].set_xticks(x,names,rotation=60,ha='right');fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        for track,names,corr in heatmaps:
            fig,ax=plt.subplots(figsize=(9,8));im=ax.imshow(corr,vmin=-1,vmax=1,cmap='coolwarm')
            ax.set_xticks(range(len(names)),names,rotation=60,ha='right');ax.set_yticks(range(len(names)),names)
            ax.set_title(f'Gene signal Spearman correlation\n{track}');fig.colorbar(im,ax=ax,label='Spearman r')
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)


if __name__=='__main__': main(snakemake)
