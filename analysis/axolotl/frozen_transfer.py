#!/usr/bin/env python3
"""Training-frozen axolotl positional-state analysis.

Reproduces the primary v0.16/v0.17 scoring logic from public processed count
matrices. All gene selection, directions, means and SDs are fit only in the
training series before projection into validation/perturbation samples.
"""
from __future__ import annotations
import argparse, itertools, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

MATURE_A = ["../HISAT2/82226_All.sort.bam", "../HISAT2/82228_All.sort.bam"]
MATURE_P = ["../HISAT2/82221_All.sort.bam", "../HISAT2/82227_All.sort.bam", "../HISAT2/82229_All.sort.bam"]
HI_A = [f"AnteriorBlastemaHI_Rep{i}" for i in range(1,5)]
HI_P = [f"PosteriorBlastemaHI_Rep{i}" for i in range(1,5)]
LOW_A = [f"AnteriorBlastemaLOW_Rep{i}" for i in range(1,5)]
LOW_P = [f"PosteriorBlastemaLOW_Rep{i}" for i in range(1,5)]
HAND2 = [f"Prrx1-Hand2_Rep{i}" for i in range(1,4)]
MCHERRY = [f"Prrx1-mCherry-A_Rep{i}" for i in range(1,3)]
ATOP = ["APtransplant_Rep2", "APtransplant_Rep3"]
ATOA = ["Aatransplant_Rep1", "Aatransplant_Rep2"]

def logcpm(df):
    a=df.to_numpy(dtype=float)
    lib=a.sum(axis=0,keepdims=True)
    lib[lib==0]=1
    return pd.DataFrame(np.log2(1e6*a/lib+1),index=df.index,columns=df.columns)

def build_axis(expr,pos,neg,k=100,expr_threshold=1.0):
    effect=expr[pos].mean(axis=1)-expr[neg].mean(axis=1)
    mean=(expr[pos].mean(axis=1)+expr[neg].mean(axis=1))/2
    diffs=np.column_stack([(expr[p]-expr[n]).to_numpy() for p in pos for n in neg])
    ev=effect.to_numpy(); valid=mean.to_numpy()>expr_threshold
    up=np.flatnonzero(valid & (diffs>0).all(axis=1))
    dn=np.flatnonzero(valid & (diffs<0).all(axis=1))
    up=up[np.argsort(ev[up])[::-1]][:k]
    dn=dn[np.argsort(ev[dn])][:k]
    idx=np.r_[up,dn]
    signs=np.r_[np.ones(len(up)), -np.ones(len(dn))]
    train=expr[pos+neg].to_numpy()[idx,:].T
    mu=train.mean(axis=0)
    sd=train.std(axis=0,ddof=1); sd[sd==0]=1
    weights=np.abs(ev[idx])
    return dict(idx=idx,signs=signs,mu=mu,sd=sd,weights=weights,effect=ev[idx])

def score_axis(axis,expr,positive,negative,weighting="sign_only"):
    cols=positive+negative
    X=expr[cols].to_numpy()[axis['idx'],:].T
    z=(X-axis['mu'])/axis['sd']
    if weighting=="sign_only":
        scores=(z*axis['signs']).mean(axis=1)
    elif weighting=="effect_weighted":
        w=axis['weights']
        scores=(z*axis['signs']*w).sum(axis=1)/w.sum()
    else: raise ValueError(weighting)
    y=np.r_[np.ones(len(positive)),np.zeros(len(negative))]
    return scores,float(roc_auc_score(y,scores))

def exact_one_sided_p(scores,n_pos):
    obs=scores[:n_pos].mean()-scores[n_pos:].mean()
    vals=[]; n=len(scores)
    for pos_idx in itertools.combinations(range(n),n_pos):
        mask=np.zeros(n,dtype=bool); mask[list(pos_idx)]=True
        vals.append(scores[mask].mean()-scores[~mask].mean())
    vals=np.asarray(vals)
    return float((vals>=obs-1e-12).mean()),len(vals)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mature',required=True)
    ap.add_argument('--regen',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args(); out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    m=pd.read_csv(args.mature,sep='\t').set_index('Geneid')
    r=pd.read_csv(args.regen).set_index('Unnamed: 0')
    common=m.index.intersection(r.index); m=m.loc[common]; r=r.loc[common]
    lm,lr=logcpm(m),logcpm(r)
    rows=[]
    mature=build_axis(lm,MATURE_P,MATURE_A)
    comps=[('blastema_HI',HI_P,HI_A),('blastema_LOW',LOW_P,LOW_A),('Hand2',HAND2,MCHERRY),('AtoP',ATOP,ATOA)]
    for w in ['sign_only','effect_weighted']:
        for name,pos,neg in comps:
            sc,auc=score_axis(mature,lr,pos,neg,w); p,npart=exact_one_sided_p(sc,len(pos))
            rows.append(dict(axis='mature',weighting=w,comparison=name,auc=auc,exact_one_sided_p=p,number_exact_partitions=npart,positive_mean=float(sc[:len(pos)].mean()),negative_mean=float(sc[len(pos):].mean()),mean_difference=float(sc[:len(pos)].mean()-sc[len(pos):].mean()),positive_scores='|'.join(map(lambda x:f'{x:.9g}',sc[:len(pos)])),negative_scores='|'.join(map(lambda x:f'{x:.9g}',sc[len(pos):]))))
    acute=build_axis(lr,HI_P,HI_A)
    comps2=[('blastema_LOW',LOW_P,LOW_A),('Hand2',HAND2,MCHERRY),('AtoP',ATOP,ATOA)]
    for w in ['sign_only','effect_weighted']:
        for name,pos,neg in comps2:
            sc,auc=score_axis(acute,lr,pos,neg,w); p,npart=exact_one_sided_p(sc,len(pos))
            rows.append(dict(axis='acute_HI',weighting=w,comparison=name,auc=auc,exact_one_sided_p=p,number_exact_partitions=npart,positive_mean=float(sc[:len(pos)].mean()),negative_mean=float(sc[len(pos):].mean()),mean_difference=float(sc[:len(pos)].mean()-sc[len(pos):].mean()),positive_scores='|'.join(map(lambda x:f'{x:.9g}',sc[:len(pos)])),negative_scores='|'.join(map(lambda x:f'{x:.9g}',sc[len(pos):]))))
    scores=pd.DataFrame(rows); scores.to_csv(out/'frozen_transfer.csv',index=False)
    # Training-frozen mature leave-one-out, primary cross-series columns only.
    loo=[]
    for omitted in MATURE_A+MATURE_P:
        p=[x for x in MATURE_P if x!=omitted]; a=[x for x in MATURE_A if x!=omitted]
        ax=build_axis(lm,p,a)
        row={'omitted_sample':omitted.split('/')[-1].split('_')[0]}
        for name,pos,neg in [('blastema_HI',HI_P,HI_A),('blastema_LOW',LOW_P,LOW_A)]:
            sc,auc=score_axis(ax,lr,pos,neg,'sign_only'); row[name+'_auc']=auc
        loo.append(row)
    loo=pd.DataFrame(loo); loo.to_csv(out/'mature_leave_one_out_training_frozen.csv',index=False)
    summary={
      'mature_HI_auc':float(scores.query("axis=='mature' and weighting=='sign_only' and comparison=='blastema_HI'").auc.iloc[0]),
      'mature_LOW_auc':float(scores.query("axis=='mature' and weighting=='sign_only' and comparison=='blastema_LOW'").auc.iloc[0]),
      'mature_Hand2_auc':float(scores.query("axis=='mature' and weighting=='sign_only' and comparison=='Hand2'").auc.iloc[0]),
      'mature_AtoP_auc':float(scores.query("axis=='mature' and weighting=='sign_only' and comparison=='AtoP'").auc.iloc[0]),
      'acute_AtoP_sign_auc':float(scores.query("axis=='acute_HI' and weighting=='sign_only' and comparison=='AtoP'").auc.iloc[0]),
      'acute_Hand2_sign_auc':float(scores.query("axis=='acute_HI' and weighting=='sign_only' and comparison=='Hand2'").auc.iloc[0]),
      'acute_Hand2_weighted_auc':float(scores.query("axis=='acute_HI' and weighting=='effect_weighted' and comparison=='Hand2'").auc.iloc[0]),
      'loo_HI_aucs':loo.blastema_HI_auc.tolist(),
      'loo_LOW_aucs':loo.blastema_LOW_auc.tolist(),
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    # Hard checks for manuscript-facing core numbers.
    assert summary['mature_HI_auc']==1.0 and summary['mature_LOW_auc']==1.0
    assert summary['mature_Hand2_auc']==1.0 and summary['mature_AtoP_auc']==0.75
    assert summary['acute_AtoP_sign_auc']==1.0
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
