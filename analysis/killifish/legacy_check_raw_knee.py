import h5py, numpy as np, pandas as pd, os, math
from pathlib import Path

files=[('/mnt/data/gse260629_h5/GSM8121285_L45793.h5',7938,'prox_rep2'),('/mnt/data/gse260629_h5/GSM8121286_L45794.h5',9060,'dist_rep2')]
rows=[]
for fn,nexp,label in files:
    with h5py.File(fn,'r') as f:
        d=f['matrix/data'][:]
        indptr=f['matrix/indptr'][:]
        # per-column sum without constructing sparse matrix
        cs=np.add.reduceat(d, indptr[:-1], dtype=np.int64)
        empty=(indptr[1:]==indptr[:-1])
        cs[empty]=0
    order=np.sort(cs)[::-1]
    for mult in [0.5,0.75,1.0,1.25,1.5,2.0,3.0]:
        k=min(len(order)-1,max(1,int(round(nexp*mult))))
        rows.append([label,nexp,mult,k,int(order[k-1])])
    q={q:float(np.quantile(order[:nexp],q)) for q in [0,.01,.05,.1,.25,.5,.75,.9,.99,1]}
    print('\n',label,'droplets',len(order),'expected',nexp)
    print('count at expected rank',int(order[nexp-1]),'next',int(order[nexp]))
    print('top-N quantiles',q)
    # largest local ratio/drop in region expected/3..expected*3 using log diff
    lo=max(10,nexp//3); hi=min(len(order)-1,nexp*3)
    arr=order[lo:hi].astype(float)+1
    drops=np.log(arr[:-1])-np.log(arr[1:])
    idx=lo+int(np.argmax(drops))
    print('largest log drop near expected at rank',idx+1,'counts',int(order[idx]),'->',int(order[idx+1]),'ratio',float((order[idx]+1)/(order[idx+1]+1)))

pd.DataFrame(rows,columns=['sample','expected_cells','rank_multiple','rank','UMI_at_rank']).to_csv('/mnt/data/raw_knee_summary.csv',index=False)
