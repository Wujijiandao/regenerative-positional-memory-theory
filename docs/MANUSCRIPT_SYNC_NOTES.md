# Manuscript / repository synchronization notes

## v1.0.0 synchronization state

The public repository and the submission manuscript have been synchronized to the **training-frozen** axolotl scoring convention.

Canonical mature leave-one-out values used by Figure 4 are:

- HI blastema AUC: 1.0 in 5/5 rebuilt signatures;
- LOW blastema AUCs: 0.75, 1.0, 0.875, 1.0, 0.875 for omitted samples 82226, 82228, 82221, 82227 and 82229, respectively.

The corresponding manuscript synchronization patch is internal project version **v0.17.1**. No theorem, dataset, biological interpretation or qualitative conclusion changed relative to v0.17.0; only the Figure 4 visualization was synchronized to the already-adopted training-frozen method.

The killifish remote-dorsal negative control is reproduced end-to-end by `analysis/killifish/remote_dorsal_control.py`. Older cell-level scripts are retained only as explicitly labelled legacy/descriptive provenance.
