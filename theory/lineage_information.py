"""Exact finite-distribution demonstrations for fixed-target lineage information.

The random variable of interest is a fixed pre-injury target Theta.

A noisy surviving mother-state S may carry less than one bit about Theta.
Copying S into arbitrarily many descendants does not increase the descendants'
joint mutual information with Theta, although the sum of marginal MIs can be
arbitrarily large.
"""

from collections import defaultdict
from itertools import product
from math import log2
from typing import Dict, Hashable, Iterable, Tuple


Prob = Dict[Tuple[Hashable, ...], float]


def entropy_from_probs(probs: Iterable[float]) -> float:
    h = 0.0
    for p in probs:
        if p > 0:
            h -= p * log2(p)
    return h


def marginal(joint: Prob, indices: Tuple[int, ...]) -> Prob:
    out = defaultdict(float)
    for state, p in joint.items():
        key = tuple(state[i] for i in indices)
        out[key] += p
    return dict(out)


def mutual_information(joint: Prob, x_idx: Tuple[int, ...], y_idx: Tuple[int, ...]) -> float:
    pxy = marginal(joint, x_idx + y_idx)
    px = marginal(joint, x_idx)
    py = marginal(joint, y_idx)

    mi = 0.0
    for xy, p in pxy.items():
        x = xy[:len(x_idx)]
        y = xy[len(x_idx):]
        if p > 0:
            mi += p * log2(p / (px[x] * py[y]))
    return mi


def noisy_seed_distribution(p_flip: float, descendants: int) -> Prob:
    """Joint state: (Theta, S, Y1, ..., Yn), with Yi=S deterministically."""
    if not (0 <= p_flip <= 1):
        raise ValueError("p_flip must be in [0,1]")
    if descendants < 1:
        raise ValueError("descendants must be >=1")

    joint = defaultdict(float)
    for theta in (0, 1):
        p_theta = 0.5
        for noise in (0, 1):
            p_n = (1 - p_flip) if noise == 0 else p_flip
            s = theta ^ noise
            ys = (s,) * descendants
            joint[(theta, s) + ys] += p_theta * p_n
    return dict(joint)


def target_seed_mi(p_flip: float) -> float:
    joint = noisy_seed_distribution(p_flip, 1)
    return mutual_information(joint, (0,), (1,))


def target_descendant_joint_mi(p_flip: float, descendants: int) -> float:
    joint = noisy_seed_distribution(p_flip, descendants)
    y_idx = tuple(range(2, 2 + descendants))
    return mutual_information(joint, (0,), y_idx)


def target_descendant_marginal_sum(p_flip: float, descendants: int) -> float:
    joint = noisy_seed_distribution(p_flip, descendants)
    return sum(
        mutual_information(joint, (0,), (2 + i,))
        for i in range(descendants)
    )


def binary_entropy(p: float) -> float:
    if p in (0.0, 1.0):
        return 0.0
    return -p * log2(p) - (1 - p) * log2(1 - p)


def theoretical_seed_mi(p_flip: float) -> float:
    return 1.0 - binary_entropy(p_flip)


def demo() -> str:
    p = 0.20
    n = 64

    seed = target_seed_mi(p)
    joint = target_descendant_joint_mi(p, n)
    marginal_sum = target_descendant_marginal_sum(p, n)

    lines = [
        "Growing-lineage fixed-target information toy",
        "--------------------------------------------",
        f"noise in surviving mother-state p={p}",
        f"I(Theta; mother state) = {seed:.6f} bits",
        f"theory 1-h2(p) = {theoretical_seed_mi(p):.6f} bits",
        f"number of daughter copies = {n}",
        f"I(Theta; all daughters jointly) = {joint:.6f} bits",
        f"sum_i I(Theta; daughter_i) = {marginal_sum:.6f} bits",
        f"joint information increased by copying = {joint > seed + 1e-10}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
