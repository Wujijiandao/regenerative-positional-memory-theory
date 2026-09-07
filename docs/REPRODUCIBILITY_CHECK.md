# Reproducibility check — release candidate 1.0.0-rc1

Automated check performed before repository upload:

- theory: 26/26 unit tests passed;
- axolotl: recomputed from GSE243137 and GSE284768 processed matrices;
- mature -> HI AUC = 1.0;
- mature -> LOW AUC = 1.0;
- mature -> Hand2 AUC = 1.0;
- mature -> A->P AUC = 0.75;
- acute HI -> A->P AUC = 1.0;
- acute HI -> Hand2 AUC = 0.50 (sign-only) and 0.333... (effect-weighted);
- training-frozen mature LOO HI AUC = 1.0 in 5/5;
- training-frozen mature LOO LOW AUCs = 0.75, 1.0, 0.875, 1.0, 0.875;
- killifish remote-dorsal 50K proximal-minus-distal score difference = 0.5004758459532846;
- killifish exact one-sided 3-vs-3 permutation p = 0.05;
- 30K/50K score Spearman rho = 0.9833333333333333;
- caudal-signature effect rho with dorsal 50K effect = 0.53822395559889;
- signature direction agreement = 0.63.

These automated checks do not replace the author's independent scientific verification required before journal submission.
