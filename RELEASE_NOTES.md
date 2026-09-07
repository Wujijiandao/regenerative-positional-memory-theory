This release freezes the reproducibility materials supporting the manuscript **“Spatial robustness and environmental side information in regenerative positional memory.”**

It includes:

- code and exhaustive tests for the connected-lesion theory and robustness–latency results;
- the total-distance-domination baseline used for the single-cell limit;
- training-frozen axolotl analyses based on GEO series **GSE243137** and **GSE284768**;
- an end-to-end African killifish remote-dorsal negative-control analysis based on **GSE260629**;
- derived tables and manuscript-facing figures;
- data manifests, sample mappings, source filenames and SHA-256 checksums;
- reproducibility instructions and environment specifications.

The biological source datasets are not redistributed. They remain available from their original public repositories and are referenced by accession and verified source metadata.

The release reproduces the principal manuscript-facing results, including:

- mature positional-state → high-input blastema: ROC AUC = 1.00;
- mature positional-state → low-input blastema: ROC AUC = 1.00;
- Hand2 perturbation projected onto the mature positional-state axis: ROC AUC = 1.00;
- anterior-to-posterior transplantation projected onto the mature positional-state axis: ROC AUC = 0.75;
- killifish remote-dorsal proximal–distal score difference ≈ 0.5005 with an exact one-sided permutation p = 0.05;
- 30K/50K remote-dorsal technical reproducibility: Spearman ρ ≈ 0.983.

This repository is intended to provide the reproducibility record associated with the manuscript submission and to be archived as a versioned DOI-bearing release through Zenodo.

### Scope and interpretation

The single-cell multi-source limit is explicitly treated as a known total-distance-domination baseline rather than a new graph-theoretic parameter. The information-theoretic bound is presented as an application of classical side-information reasoning.

The empirical analyses test the distinction between persistent local positional state and target-correlated environmental information; they do not constitute a direct experimental validation of the exact connected-lesion geometric law.
