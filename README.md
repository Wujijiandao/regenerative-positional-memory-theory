# Spatial robustness and environmental side information in regenerative positional memory

Reproducibility code and derived results supporting the manuscript:

**Spatial robustness and environmental side information in regenerative positional memory**  
Yuzhan Zhang — Independent researcher

This repository accompanies a quantitative framework for regenerative positional information under tissue injury. The manuscript separates three questions that are often conflated: (i) whether a connected lesion can erase target-associated state, (ii) how far surviving state must propagate to a reconstituting region, and (iii) how much target-correlated information can be supplied by the post-injury environment.

## What is in this repository

- `theory/` — graph-theoretic and information-theoretic model code plus exhaustive finite checks;
- `analysis/axolotl/` — training-frozen reanalysis of GSE243137 and GSE284768;
- `analysis/killifish/` — caudal-signature / remote-dorsal negative-control analysis for GSE260629;
- `results/` — manuscript-facing derived tables and validation summaries;
- `figures/` — frozen main and supplementary figures;
- `data_manifest/` — source accessions, filenames, SHA-256 hashes and sample mappings;
- `scripts/` — data verification and reproduction helpers;
- `.github/workflows/` — lightweight theory test workflow that does not require third-party biological data.

## Main scope and caveats

The repository does **not** claim that total-distance domination, unreliable facility location or side-information theory are new mathematical ideas. The single-cell backup limit is treated as a known total-distance-domination baseline. The application-specific theory concerns connected tissue lesions, target-state erasure and post-lesion access geometry.

The public biological datasets test the proposed **information-source decomposition**, not the exact geometric law itself. Axolotl transfer is a same-study cross-series analysis, not an independent-laboratory replication. Killifish cell-level analyses are descriptive where independent biological replication is unavailable.

## Source data

Third-party raw/processed source data are not redistributed by this repository. Obtain them from the original repositories and place them under `data/` as described in `data_manifest/README.md`.

Primary accessions:

- **GSE260629** — African killifish fin regeneration;
- **GSE243137** — uninjured anterior/posterior axolotl dermal connective tissue;
- **GSE284768** — axolotl blastema, Hand2 perturbation and transplantation experiments.

The repository records expected filenames and SHA-256 hashes so that the exact source snapshots used in the analysis can be checked.

## Quick start

Create an environment:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\\Scripts\\Activate.ps1

pip install -r requirements.txt
```

Run the theory checks:

```bash
cd theory
python test_models.py
```

Verify downloaded data:

```bash
python scripts/verify_data.py --data-dir data
```

Run the training-frozen axolotl analysis after placing the two processed matrices in `data/axolotl/`:

```bash
python analysis/axolotl/frozen_transfer.py \
  --mature data/axolotl/GSE243137_MatrixCounts_allsamples.txt.gz \
  --regen data/axolotl/GSE284768_MatrixCounts.csv.gz \
  --out results/axolotl/recomputed
```

Run the killifish remote-dorsal negative control after downloading `GSE260629_RAW.tar`:

```bash
python analysis/killifish/remote_dorsal_control.py \
  --geo-tar data/killifish/GSE260629_RAW.tar \
  --out results/killifish/recomputed
```

## Reproducibility status

This repository is prepared as a **submission/reviewer reproducibility release candidate**. Before a public `v1.0.0` archival release, the author should independently rerun the complete workflow, verify manuscript numbers and figures, and finalize the data/code citation.

See `docs/HUMAN_VERIFICATION.md` and `docs/MANUSCRIPT_SYNC_NOTES.md`.

## License

Original code and repository-authored documentation are released under the MIT License. Third-party source datasets and supplementary files are **not** included and remain subject to their original terms.

## Citation

A Zenodo DOI will be inserted here after the archival `v1.0.0` release. Until then, use the metadata in `CITATION.cff`.
