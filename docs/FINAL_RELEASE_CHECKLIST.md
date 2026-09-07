# Final release checklist — v1.0.0

Automated repository QA:

- [x] transient `__pycache__` / `.pyc` files removed;
- [x] theory unit tests pass;
- [x] axolotl analysis recomputed from frozen processed matrices;
- [x] killifish remote-dorsal analysis recomputed from GSE260629_RAW.tar;
- [x] Figure 4 synchronized to training-frozen leave-one-out values;
- [x] source-data filenames and SHA-256 hashes recorded;
- [x] `CITATION.cff` version set to 1.0.0;
- [x] MIT licence applies only to repository-authored code/documentation;
- [x] third-party source data are not redistributed;
- [x] GitHub release notes and Zenodo metadata are frozen.

Release operations:

1. push the frozen repository tree to `main`;
2. create GitHub tag/release `v1.0.0`;
3. archive the GitHub release in Zenodo;
4. add the minted Zenodo DOI back to `README.md` and `CITATION.cff` in a DOI-metadata follow-up commit if desired;
5. use the DOI in the journal Data Accessibility statement.

Independent author verification remains required before journal submission; it is not represented by this automated checklist.
