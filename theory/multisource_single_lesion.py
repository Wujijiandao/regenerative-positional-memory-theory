"""Exact multi-source theorem for a one-dimensional tissue under one-cell injury.

For m >= 2 identical target-memory carriers on a path of n cells:

    R_m*(n,1)
      = max(1, ceil((n-floor(m/2)) / (m+floor(m/2))))

The code provides:
- exact formula;
- brute-force validation for small systems;
- a constructive modular placement.
"""

from itertools import combinations
from math import ceil
from typing import Tuple, List


def worst_radius_single_cell(n: int, sources: Tuple[int, ...]) -> int:
    S = set(sources)
    if len(S) < 2:
        return float("inf")
    worst = 0
    for x in range(1, n + 1):
        surviving = S - {x}
        d = min(abs(x - s) for s in surviving)
        worst = max(worst, d)
    return worst


def theorem_radius(n: int, m: int) -> int:
    if not (2 <= m <= n):
        raise ValueError("Require 2 <= m <= n")
    c = m // 2
    return max(1, ceil((n - c) / (m + c)))


def ordinary_no_failure_radius(n: int, m: int) -> int:
    if not (1 <= m <= n):
        raise ValueError("Require 1 <= m <= n")
    return max(0, ceil((n - m) / (2 * m)))


def brute_force_optimum(n: int, m: int):
    best = float("inf")
    placements = []
    for S in combinations(range(1, n + 1), m):
        r = worst_radius_single_cell(n, S)
        if r < best:
            best = r
            placements = [S]
        elif r == best:
            placements.append(S)
    return best, placements


def _decompose(total: int, lows: List[int], highs: List[int]) -> List[int]:
    """Choose integers x_i in [low_i, high_i] summing to total."""
    if total < sum(lows) or total > sum(highs):
        raise ValueError("total outside feasible range")
    out = lows[:]
    rem = total - sum(out)
    for i in range(len(out)):
        add = min(rem, highs[i] - lows[i])
        out[i] += add
        rem -= add
    if rem:
        raise RuntimeError("decomposition failed")
    return out


def _place_module(start: int, length: int, k: int, r: int) -> List[int]:
    """Place k=2 or 3 sources in a block [start, start+length-1]."""
    if k not in (2, 3):
        raise ValueError("module size must be 2 or 3")
    target = length - 1

    # Gaps: left endpoint, internal k-1 gaps, right endpoint.
    # Internal gaps start at 1; all capacities are <= r.
    gaps = [0] + [1] * (k - 1) + [0]
    caps = [r] + [r] * (k - 1) + [r]
    current = sum(gaps)
    rem = target - current
    for i in range(len(gaps)):
        room = caps[i] - gaps[i]
        add = min(rem, room)
        gaps[i] += add
        rem -= add
    if rem != 0:
        raise RuntimeError("could not realize module length")

    pos = start + gaps[0]
    sources = [pos]
    for d in gaps[1:-1]:
        pos += d
        sources.append(pos)
    return sources


def construct_optimal(n: int, m: int) -> Tuple[int, ...]:
    r = theorem_radius(n, m)
    c = m // 2

    if m % 2 == 0:
        sizes = [2] * c
    else:
        sizes = [2] * (c - 1) + [3]

    lows = sizes[:]
    highs = [(k + 1) * r + 1 for k in sizes]
    lengths = _decompose(n, lows, highs)

    sources = []
    start = 1
    for length, k in zip(lengths, sizes):
        sources.extend(_place_module(start, length, k, r))
        start += length

    result = tuple(sorted(sources))
    if len(result) != m:
        raise RuntimeError("wrong number of sources")
    if worst_radius_single_cell(n, result) > r:
        raise RuntimeError("construction failed")
    return result


def demo() -> str:
    lines = [
        "Exact multi-source single-cell-injury theorem",
        "--------------------------------------------",
    ]
    for n, m in [(40, 4), (60, 6), (80, 9), (120, 12)]:
        r = theorem_radius(n, m)
        S = construct_optimal(n, m)
        r0 = ordinary_no_failure_radius(n, m)
        lines += [
            f"n={n}, m={m}:",
            f"  robust optimum R*={r}",
            f"  no-failure optimum R0={r0}",
            f"  constructed carriers={S}",
            f"  verified construction radius={worst_radius_single_cell(n,S)}",
        ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
