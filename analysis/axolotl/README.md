# Axolotl analysis

Primary source series: GSE243137 and GSE284768 (same published study).

`frozen_transfer.py` fits gene selection, direction, mean and standard deviation only in the training series and then projects validation/perturbation samples without refitting normalization parameters in the validation cohort.

Expected core results under sign-only frozen scoring:
- mature axis -> blastema HI: AUC 1.0;
- mature axis -> blastema LOW: AUC 1.0;
- mature axis -> Hand2 vs mCherry: AUC 1.0;
- mature axis -> A->P vs A->A: AUC 0.75;
- acute HI axis -> A->P vs A->A: AUC 1.0;
- acute HI axis -> Hand2: AUC 0.50 sign-only and 0.333 effect-weighted.

The script also recomputes the mature-sample leave-one-out transfer using training-frozen scaling.
