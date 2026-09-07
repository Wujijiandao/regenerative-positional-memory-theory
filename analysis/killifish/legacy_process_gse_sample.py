import h5py, numpy as np, scipy.sparse as sp, pandas as pd, sys, json
from pathlib import Path

MITO = {f'ENSNFUG000150000{i:02d}' for i in range(2,40)}

def load_h5(fn, expected=None):
    with h5py.File(fn,'r') as f:
        g=f['matrix']
        data=g['data'][:]
        indices=g['indices'][:]
        indptr=g['indptr'][:]
        shape=tuple(g['shape'][:])
        ids=np.array(g['features/id']).astype('U')
        names=np.array(g['features/name']).astype('U')
        ft=np.array(g['features/feature_type']).astype('U') if 'feature_type' in g['features'] else np.array(['Gene Expression']*shape[0])
        barcodes=np.array(g['barcodes']).astype('U')
    X=sp.csc_matrix((data,indices,indptr),shape=shape)
    gex=np.flatnonzero(ft=='Gene Expression')
    X=X[gex,:]
    ids=ids[gex]; names=names[gex]
    if expected is not None and X.shape[1] > expected*10:
        totals=np.asarray(X.sum(axis=0)).ravel()
        keep=np.argpartition(totals,-expected)[-expected:]
        keep=np.sort(keep)
        X=X[:,keep]
        barcodes=barcodes[keep]
    return X,ids,names,barcodes

def qc(X,ids):
    total=np.asarray(X.sum(axis=0)).ravel().astype(float)
    nfeat=np.diff(X.indptr).astype(float)
    mitorows=np.flatnonzero(np.isin(ids,list(MITO)))
    mt=np.asarray(X[mitorows,:].sum(axis=0)).ravel().astype(float)
    pct=np.divide(100*mt,total,out=np.zeros_like(mt),where=total>0)
    mp=pct[pct<40]
    mito_cut=np.median(mp)+2*np.std(mp,ddof=1)
    nf=nfeat[pct<mito_cut]
    med=np.median(nf); sd=np.std(nf,ddof=1)
    low=med-1.2*sd; high=med+4*sd
    keep=(nfeat>low)&(nfeat<high)&(pct<mito_cut)
    return keep, dict(input_cells=int(X.shape[1]),retained_cells=int(keep.sum()),mito_cut=float(mito_cut),feature_low=float(low),feature_high=float(high),median_umi=float(np.median(total[keep])),median_features=float(np.median(nfeat[keep])),median_mt=float(np.median(pct[keep])))

fn=sys.argv[1]; out=sys.argv[2]; expected=None if sys.argv[3]=='None' else int(sys.argv[3])
X,ids,names,barcodes=load_h5(fn,expected)
keep,stats=qc(X,ids)
Xq=X[:,keep].tocsr()
barq=barcodes[keep]
sp.save_npz(out+'.npz',Xq)
np.save(out+'_ids.npy',ids)
np.save(out+'_names.npy',names)
np.save(out+'_barcodes.npy',barq)
Path(out+'_qc.json').write_text(json.dumps(stats,indent=2))
print(json.dumps(stats))
