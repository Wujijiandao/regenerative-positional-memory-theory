"""Exact two-source robustness–latency tradeoff on a path.

Tissue: vertices 1..n.
Two target-memory sources: a < b.
Allowed injury: any connected interval of length <= L.

Survival requires b-a >= L.

Worst access radius:
    max over allowed lesions D
    max over x in D
    distance(x, surviving sources).

For n >= 3L+1, the exact optimum over two source placements is

    R* = max(2L-1, ceil((n+L-2)/3)).

This module brute-forces the injury model and verifies the theorem.
"""

from itertools import combinations
from math import ceil
from typing import Tuple, List


def exact_access_formula(n: int, a: int, b: int, L: int) -> int:
    if not (1 <= a < b <= n):
        raise ValueError("Require 1 <= a < b <= n")
    if L < 1:
        raise ValueError("L must be positive")
    if b - a < L:
        return float("inf")

    return max(
        a - 1,
        n - b,
        b - max(1, a - L + 1),
        min(n, b + L - 1) - a,
    )


def brute_force_access_radius(n: int, a: int, b: int, L: int) -> int:
    sources = {a, b}
    worst = 0
    for left in range(1, n + 1):
        for right in range(left, min(n, left + L - 1) + 1):
            D = set(range(left, right + 1))
            surviving = sources - D
            if not surviving:
                return float("inf")
            local = max(
                min(abs(x - s) for s in surviving)
                for x in D
            )
            worst = max(worst, local)
    return worst


def theorem_optimum(n: int, L: int) -> int:
    if n < 3 * L + 1:
        raise ValueError("Closed-form theorem assumes n >= 3L+1")
    return max(
        2 * L - 1,
        ceil((n + L - 2) / 3),
    )


def brute_force_optimum(n: int, L: int):
    best = float("inf")
    placements = []
    for a in range(1, n):
        for b in range(a + 1, n + 1):
            if b - a < L:
                continue
            r = brute_force_access_radius(n, a, b, L)
            if r < best:
                best = r
                placements = [(a, b)]
            elif r == best:
                placements.append((a, b))
    return best, placements


def construction_window(n: int, L: int):
    R = theorem_optimum(n, L)
    low_s = max(L, n - 1 - 2 * R)
    high_s = R - L + 1
    return R, low_s, high_s


def construct_optimal_pair(n: int, L: int) -> Tuple[int, int]:
    R, low_s, high_s = construction_window(n, L)
    if low_s > high_s:
        raise RuntimeError("construction window unexpectedly empty")

    s = low_s
    # Need a-1 <= R and n-(a+s) <= R.
    low_a = max(1, n - s - R)
    high_a = min(R + 1, n - s)
    if low_a > high_a:
        raise RuntimeError("no centered placement exists")
    a = low_a
    b = a + s

    if exact_access_formula(n, a, b, L) > R:
        raise RuntimeError("construction failed")
    return a, b


def demo() -> str:
    examples = [(31, 4), (40, 5), (61, 8)]
    lines = [
        "Two-source positional-memory robustness–latency theorem",
        "------------------------------------------------------",
    ]
    for n, L in examples:
        R = theorem_optimum(n, L)
        a, b = construct_optimal_pair(n, L)
        brute, _ = brute_force_optimum(n, L)
        lines += [
            f"n={n}, L={L}:",
            f"  theorem optimum R* = {R}",
            f"  constructed sources = ({a}, {b})",
            f"  connected-lesion distance = {b-a+1}",
            f"  exact access radius = {exact_access_formula(n,a,b,L)}",
            f"  brute-force optimum = {brute}",
            f"  theorem sharp = {R == brute}",
        ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())


def fixed_separation_optimum(n: int, L: int, s: int) -> int:
    """Exact fixed-separation optimum in the interior regime L <= s <= n-2L+1."""
    if not (1 <= L and L <= s <= n - 2 * L + 1):
        raise ValueError("require 1 <= L <= s <= n-2L+1")
    return max(s + L - 1, ceil((n - s - 1) / 2))

def brute_force_fixed_separation(n: int, L: int, s: int) -> int:
    if not (1 <= L <= s < n):
        raise ValueError("invalid parameters")
    vals=[]
    for a in range(1,n-s+1):
        vals.append(brute_force_access_radius(n,a,a+s,L))
    return min(vals)
