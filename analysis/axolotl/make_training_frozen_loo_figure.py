#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser()
ap.add_argument('--input',required=True)
ap.add_argument('--output',required=True)
args=ap.parse_args()
df=pd.read_csv(args.input)
mat=df[['blastema_HI_auc','blastema_LOW_auc']].to_numpy()
fig=plt.figure(figsize=(6.8,5.2))
ax=fig.add_subplot(111)
im=ax.imshow(mat,vmin=0,vmax=1,aspect='auto')
fig.colorbar(im,ax=ax,label='ROC AUC')
ax.set_xticks([0,1],['Blastema HI','Blastema LOW'])
ax.set_yticks(range(len(df)),df['omitted_sample'].astype(str))
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j,i,f'{mat[i,j]:.3g}',ha='center',va='center')
ax.set_title('Training-frozen leave-one-mature-sample-out transfer')
ax.set_ylabel('Omitted mature sample')
fig.tight_layout()
Path(args.output).parent.mkdir(parents=True,exist_ok=True)
fig.savefig(args.output,dpi=300)
plt.close(fig)
