# Killifish analysis

The manuscript-facing inferential control is `remote_dorsal_control.py`, which reconstructs the 200-gene caudal signature from CellPlex proximal/distal libraries and projects it into 3 proximal-fish and 3 distal-fish uninjured dorsal samples at both sequencing depths. The exact 3-vs-3 permutation test uses fish as biological units.

The `legacy_*` scripts preserve the earlier cell-level exploratory processing. They are retained for provenance; cell-level AUCs are descriptive and are not treated as independent biological-replicate inference.
