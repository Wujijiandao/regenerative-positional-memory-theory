# Source data manifest

This repository intentionally does not redistribute third-party biological source data. Download the named files from GEO / the original publication supplementary materials, place them at the `target` paths in `source_data.json`, and run:

```bash
python scripts/verify_data.py --data-dir data
```

The hashes identify the exact snapshots used for the archived analysis.

Axolotl sample mapping: GSE243137 anterior = 82226, 82228; posterior = 82221, 82227, 82229. GSE284768 contains 4+4 HI blastema, 4+4 LOW blastema, Hand2 3 vs mCherry 2, and A->P 2 vs A->A 2.
