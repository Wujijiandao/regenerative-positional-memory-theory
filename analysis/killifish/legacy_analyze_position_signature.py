import numpy as np, scipy.sparse as sp, pandas as pd, json, math
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score
from pathlib import Path

base=Path('/mnt/data/gse_processed')
samples=['prox_rep1','dist_rep1','prox_rep2','dist_rep2']
X={s:sp.load_npz(base/f'{s}.npz').T.tocsr() for s in samples}  # cells x genes
ids=np.load(base/'prox_rep1_ids.npy',allow_pickle=False)
names=np.load(base/'prox_rep1_names.npy',allow_pickle=False)
for s in samples:
    ids2=np.load(base/f'{s}_ids.npy',allow_pickle=False)
    assert np.array_equal(ids,ids2)

# pseudobulk CPM log2
pb={}
detect={}
for s in samples:
    sums=np.asarray(X[s].sum(axis=0)).ravel().astype(float)
    pb[s]=np.log2(1e6*sums/sums.sum()+1)
    detect[s]=np.asarray((X[s]>0).mean(axis=0)).ravel()

lfc1=pb['prox_rep1']-pb['dist_rep1']
lfc2=pb['prox_rep2']-pb['dist_rep2']
expr=(pb['prox_rep1']+pb['dist_rep1']+pb['prox_rep2']+pb['dist_rep2'])/4
mask=expr>1
rho=spearmanr(lfc1[mask],lfc2[mask]).statistic
pear=pearsonr(lfc1[mask],lfc2[mask]).statistic
sign=((lfc1[mask]>0)==(lfc2[mask]>0)).mean()

# top discovery genes from rep1 only, prevalence > 5% in either rep1 condition
prev=np.maximum(detect['prox_rep1'],detect['dist_rep1'])
elig=prev>=0.05
idx=np.flatnonzero(elig)
ord_idx=idx[np.argsort(np.abs(lfc1[idx]))[::-1]]
# top 200 excluding unnamed Ensembl-only? retain all genes for model
sel=ord_idx[:200]
# balance 100 each direction if possible
prox_up=idx[lfc1[idx]>0]
prox_up=prox_up[np.argsort(lfc1[prox_up])[::-1]][:100]
dist_up=idx[lfc1[idx]<0]
dist_up=dist_up[np.argsort(lfc1[dist_up])][:100]
sel2=np.concatenate([prox_up,dist_up])

# per-platform signed z-score module: selected from rep1 only.
def scores(sample_a,sample_b):
    A=X[sample_a][:,sel2].astype(float)
    B=X[sample_b][:,sel2].astype(float)
    # library normalize log1p
    def norm(M):
        lib=np.asarray(M.sum(axis=1)).ravel()
        scale=np.divide(1e4,lib,out=np.zeros_like(lib),where=lib>0)
        M=sp.diags(scale)@M
        M.data=np.log1p(M.data)
        return M.toarray()
    Za=norm(A); Zb=norm(B)
    Z=np.vstack([Za,Zb])
    mu=Z.mean(axis=0); sd=Z.std(axis=0,ddof=1); sd[sd==0]=1
    Z=(Z-mu)/sd
    signs=np.concatenate([np.ones(len(prox_up)),-np.ones(len(dist_up))])
    score=(Z*signs).mean(axis=1)
    y=np.concatenate([np.ones(Za.shape[0]),np.zeros(Zb.shape[0])])
    return score,y
score1,y1=scores('prox_rep1','dist_rep1')
score2,y2=scores('prox_rep2','dist_rep2')
auc1=roc_auc_score(y1,score1)
auc2=roc_auc_score(y2,score2)

# cl16/bsEp2 author marker panel explicit from Fig4
basal_ids='''ENSNFUG00015017520 ENSNFUG00015022339 ENSNFUG00015024630 ENSNFUG00015016677 ENSNFUG00015009305 ENSNFUG00015003801 ENSNFUG00015011812 ENSNFUG00015016932 ENSNFUG00015000487 ENSNFUG00015000155 ENSNFUG00015005925 ENSNFUG00015016769 ENSNFUG00015005723 ENSNFUG00015003666 ENSNFUG00015012831 ENSNFUG00015000605 ENSNFUG00015009852 ENSNFUG00015003756 ENSNFUG00015004891 ENSNFUG00015016617 ENSNFUG00015011642 ENSNFUG00015024661 ENSNFUG00015001737 ENSNFUG00015014821 ENSNFUG00015019847 ENSNFUG00015015460 ENSNFUG00015019843 ENSNFUG00015009008 ENSNFUG00015024765 ENSNFUG00015022197 ENSNFUG00015007953 ENSNFUG00015013348 ENSNFUG00015007128 ENSNFUG00015024272 ENSNFUG00015021026 ENSNFUG00015025239 ENSNFUG00015020795 ENSNFUG00015015632 ENSNFUG00015014004 ENSNFUG00015002765 ENSNFUG00015013241 ENSNFUG00015004781 ENSNFUG00015008743 ENSNFUG00015017894 ENSNFUG00015024214 ENSNFUG00015001959 ENSNFUG00015021014 ENSNFUG00015022574 ENSNFUG00015009311 ENSNFUG00015005307 ENSNFUG00015012364 ENSNFUG00015023394 ENSNFUG00015023645 ENSNFUG00015022853 ENSNFUG00015013420 ENSNFUG00015010520 ENSNFUG00015011841 ENSNFUG00015013235 ENSNFUG00015002521 ENSNFUG00015024275 ENSNFUG00015008390 ENSNFUG00015011331 ENSNFUG00015004686 ENSNFUG00015024280 ENSNFUG00015025337 ENSNFUG00015006170 ENSNFUG00015006483 ENSNFUG00015025057 ENSNFUG00015009398 ENSNFUG00015024059 ENSNFUG00015025422 ENSNFUG00015022531 ENSNFUG00015006050 ENSNFUG00015008543 ENSNFUG00015012790 ENSNFUG00015021263 ENSNFUG00015013573 ENSNFUG00015009372 ENSNFUG00015010147 ENSNFUG00015005009 ENSNFUG00015018488 ENSNFUG00015000137 ENSNFUG00015022940 ENSNFUG00015023198 ENSNFUG00015016456 ENSNFUG00015017011 ENSNFUG00015020733 ENSNFUG00015025083 ENSNFUG00015010076 ENSNFUG00015023946 ENSNFUG00015005640 ENSNFUG00015017381 ENSNFUG00015023090'''.split()
idmap={x:i for i,x in enumerate(ids)}
basal_idx=np.array([idmap[x] for x in basal_ids if x in idmap],dtype=int)

def module_mean(s,idx):
    M=X[s][:,idx].astype(float)
    lib=np.asarray(X[s].sum(axis=1)).ravel().astype(float)
    scale=np.divide(1e4,lib,out=np.zeros_like(lib,dtype=float),where=lib>0)
    M=sp.diags(scale)@M; M.data=np.log1p(M.data)
    return np.asarray(M.mean(axis=1)).ravel()
basal={s:module_mean(s,basal_idx) for s in samples}

# find sqstm1 features
sq=np.flatnonzero(np.char.find(np.char.lower(names.astype(str)),'sqstm1')>=0)

# save gene table top consistent
cons=((lfc1>0)&(lfc2>0))|((lfc1<0)&(lfc2<0))
scoremag=np.minimum(np.abs(lfc1),np.abs(lfc2))
cidx=np.flatnonzero(cons & elig)
cidx=cidx[np.argsort(scoremag[cidx])[::-1]]
tab=pd.DataFrame({
    'gene_id':ids[cidx[:300]],'gene_name':names[cidx[:300]],
    'log2FC_rep1_prox_vs_dist':lfc1[cidx[:300]],
    'log2FC_rep2_prox_vs_dist':lfc2[cidx[:300]],
    'min_abs_log2FC':scoremag[cidx[:300]],
    'detect_prox_rep1':detect['prox_rep1'][cidx[:300]],
    'detect_dist_rep1':detect['dist_rep1'][cidx[:300]],
})
tab.to_csv('/mnt/data/position_consistent_genes.csv',index=False)

# scores and basal to csv for plots
pd.DataFrame({'score':np.r_[score1,score2],
              'condition':np.r_[np.where(y1==1,'prox','dist'),np.where(y2==1,'prox','dist')],
              'replicate':np.r_[np.repeat('rep1',len(y1)),np.repeat('rep2',len(y2))]}).to_csv('/mnt/data/position_scores.csv',index=False)
rows=[]
for s in samples:
    cond='prox' if s.startswith('prox') else 'dist'; rep='rep1' if s.endswith('rep1') else 'rep2'
    for val in basal[s]: rows.append((s,cond,rep,val))
pd.DataFrame(rows,columns=['sample','condition','replicate','bsEp2_marker_score']).to_csv('/mnt/data/bsep2_scores.csv',index=False)

summary={
 'genes_compared':int(mask.sum()),
 'pseudobulk_spearman':float(rho),
 'pseudobulk_pearson':float(pear),
 'pseudobulk_sign_agreement':float(sign),
 'module_genes_prox':int(len(prox_up)),
 'module_genes_dist':int(len(dist_up)),
 'position_score_auc_rep1_discovery':float(auc1),
 'position_score_auc_rep2_validation':float(auc2),
 'author_bsEp2_marker_genes_present':int(len(basal_idx)),
 'sqstm1_matches':[{'id':ids[i],'name':names[i]} for i in sq],
 'bsEp2_score_means':{s:float(basal[s].mean()) for s in samples},
 'bsEp2_score_medians':{s:float(np.median(basal[s])) for s in samples},
}
Path('/mnt/data/position_analysis_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print('\nTop consistent genes:')
print(tab.head(20).to_string(index=False))
