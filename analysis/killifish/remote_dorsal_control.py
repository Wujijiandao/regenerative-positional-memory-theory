#!/usr/bin/env python3
"""Recompute the GSE260629 remote-dorsal negative control from the GEO RAW tar."""
from __future__ import annotations
import argparse, itertools, json, tarfile, tempfile
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import spearmanr

MITO={f'ENSNFUG000150000{i:02d}' for i in range(2,40)}
CAUDAL={'prox':'GSM8121287_L45975.h5','dist':'GSM8121288_L45976.h5'}
DORSAL=[
 ('unc',1,'30K','GSM8121289_L46866.h5'),('unc',1,'50K','GSM8121290_L46866-1.h5'),
 ('prox',1,'30K','GSM8121291_L46867.h5'),('prox',1,'50K','GSM8121292_L46867-1.h5'),
 ('dist',1,'30K','GSM8121293_L46868.h5'),('dist',1,'50K','GSM8121294_L46868-1.h5'),
 ('unc',2,'30K','GSM8121295_L46869.h5'),('unc',2,'50K','GSM8121296_L46869-1.h5'),
 ('prox',2,'30K','GSM8121297_L46870.h5'),('prox',2,'50K','GSM8121298_L46870-1.h5'),
 ('dist',2,'30K','GSM8121299_L46871.h5'),('dist',2,'50K','GSM8121300_L46871-1.h5'),
 ('unc',3,'30K','GSM8121301_L46872.h5'),('unc',3,'50K','GSM8121302_L46872-1.h5'),
 ('prox',3,'30K','GSM8121303_L46873.h5'),('prox',3,'50K','GSM8121304_L46873-1.h5'),
 ('dist',3,'30K','GSM8121305_L46874.h5'),('dist',3,'50K','GSM8121306_L46874-1.h5')]

def load_h5(fn):
    with h5py.File(fn,'r') as f:
        g=f['matrix']; X=sp.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
        ids=np.array(g['features/id']).astype('U'); names=np.array(g['features/name']).astype('U'); ft=np.array(g['features/feature_type']).astype('U')
    gex=np.flatnonzero(ft=='Gene Expression')
    return X[gex,:],ids[gex],names[gex]

def qc(X,ids):
    total=np.asarray(X.sum(axis=0)).ravel().astype(float)
    nfeat=np.diff(X.indptr).astype(float)
    mtidx=np.flatnonzero(np.isin(ids,list(MITO)))
    mt=np.asarray(X[mtidx,:].sum(axis=0)).ravel().astype(float)
    pct=np.divide(100*mt,total,out=np.zeros_like(mt,dtype=float),where=total>0)
    mp=pct[pct<40]
    mito_cut=np.median(mp)+2*np.std(mp,ddof=1)
    nf=nfeat[pct<mito_cut]
    med=np.median(nf); sd=np.std(nf,ddof=1)
    low=med-1.2*sd; high=med+4*sd
    keep=(nfeat>low)&(nfeat<high)&(pct<mito_cut)
    return X[:,keep],dict(input_cells=int(X.shape[1]),retained_cells=int(keep.sum()),mito_cut=float(mito_cut),feature_low=float(low),feature_high=float(high))

def pb_logcpm(X):
    s=np.asarray(X.sum(axis=1)).ravel().astype(float)
    return np.log2(1e6*s/s.sum()+1)

def exact_p3v3(prox,dist):
    scores=np.r_[prox,dist]; obs=np.mean(prox)-np.mean(dist); vals=[]
    for comb in itertools.combinations(range(6),3):
        m=np.zeros(6,bool); m[list(comb)]=True; vals.append(scores[m].mean()-scores[~m].mean())
    return float((np.asarray(vals)>=obs-1e-12).mean())

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--geo-tar',required=True); ap.add_argument('--out',required=True); args=ap.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    needed=set(CAUDAL.values())|{x[3] for x in DORSAL}
    with tempfile.TemporaryDirectory() as td:
        td=Path(td)
        with tarfile.open(args.geo_tar,'r') as tf:
            members={m.name:m for m in tf.getmembers() if m.name in needed}
            missing=needed-set(members)
            if missing: raise RuntimeError(f'Missing H5 members: {sorted(missing)}')
            for name in needed: tf.extract(members[name],td,filter='data')
        Xp,ids,names=load_h5(td/CAUDAL['prox']); Xp,qcp=qc(Xp,ids)
        Xd,ids2,names2=load_h5(td/CAUDAL['dist']); Xd,qcd=qc(Xd,ids2)
        if not np.array_equal(ids,ids2): raise RuntimeError('Feature mismatch')
        pb_prox,pb_dist=pb_logcpm(Xp),pb_logcpm(Xd); caudal_effect=pb_prox-pb_dist
        detp=np.asarray((Xp>0).mean(axis=1)).ravel(); detd=np.asarray((Xd>0).mean(axis=1)).ravel(); elig=np.maximum(detp,detd)>=0.05
        cand=np.flatnonzero(elig)
        up=cand[caudal_effect[cand]>0]; up=up[np.argsort(caudal_effect[up])[::-1]][:100]
        dn=cand[caudal_effect[cand]<0]; dn=dn[np.argsort(caudal_effect[dn])][:100]
        sel=np.r_[up,dn]; signs=np.r_[np.ones(len(up)),-np.ones(len(dn))]
        rows=[]; qcrows=[]; pbs={}
        for cond,rep,depth,fn in DORSAL:
            X,idsx,_=load_h5(td/fn); X,q=qc(X,idsx); pb=pb_logcpm(X); pbs[(cond,rep,depth)]=pb
            qcrows.append(dict(condition=cond,replicate=f'rep{rep}',depth=depth,file=fn,**q))
        for depth in ['30K','50K']:
            keys=[k for k in pbs if k[2]==depth]
            arr=np.vstack([pbs[k][sel] for k in keys]); mu=arr.mean(0); sd=arr.std(0,ddof=1); sd[sd==0]=1
            sc=((arr-mu)/sd*signs).mean(1)
            for k,s in zip(keys,sc): rows.append(dict(condition=k[0],replicate=f'rep{k[1]}',depth=k[2],caudal_signature_score=float(s)))
        score=pd.DataFrame(rows); score.to_csv(out/'dorsal_caudal_signature_scores.csv',index=False)
        pd.DataFrame(qcrows).to_csv(out/'dorsal_qc.csv',index=False)
        tests={}
        for depth in ['30K','50K']:
            z=score[score.depth==depth]; prox=z[z.condition=='prox'].caudal_signature_score.to_numpy(); dist=z[z.condition=='dist'].caudal_signature_score.to_numpy()
            tests[depth]=dict(prox_mean=float(prox.mean()),dist_mean=float(dist.mean()),mean_difference=float(prox.mean()-dist.mean()),exact_one_sided_p=exact_p3v3(prox,dist))
        a=score[score.depth=='30K'].sort_values(['condition','replicate']).caudal_signature_score.to_numpy(); b=score[score.depth=='50K'].sort_values(['condition','replicate']).caudal_signature_score.to_numpy()
        tech=float(spearmanr(a,b).statistic)
        prox50=np.mean([pbs[('prox',i,'50K')] for i in [1,2,3]],axis=0); dist50=np.mean([pbs[('dist',i,'50K')] for i in [1,2,3]],axis=0); dorsal_effect=prox50-dist50
        sig_rho=float(spearmanr(caudal_effect[sel],dorsal_effect[sel]).statistic); sig_dir=float(np.mean(np.sign(caudal_effect[sel])==np.sign(dorsal_effect[sel])))
        gene=pd.DataFrame({'gene_id':ids,'gene_name':names,'caudal_rep1_log2FC_prox_minus_dist':caudal_effect,'dorsal_50K_groupmean_log2FC_proxFish_minus_distFish':dorsal_effect,'in_caudal_200gene_signature':np.isin(np.arange(len(ids)),sel)})
        gene.to_csv(out/'caudal_vs_dorsal_gene_effects.csv',index=False)
        summary={'biological_units':'3 independent fish per caudal condition in uncut dorsal fin','technical_depths':['30K','50K'],'primary_depth':'50K','caudal_signature_genes':int(len(sel)),'score_tests':tests,'technical_30K_50K_spearman':tech,'signature_cauda_vs_dorsal_effect_spearman_50K':sig_rho,'signature_direction_agreement_cauda_vs_dorsal_50K':sig_dir}
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        # Exact v0.11 score-level checks.
        assert abs(tests['50K']['mean_difference']-0.5004758459532846)<1e-9
        assert tests['50K']['exact_one_sided_p']==0.05
        assert abs(tech-0.9833333333333333)<1e-12
        print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
