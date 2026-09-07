"""Geometric robustness amplification without target-information gain.

A target bit is copied into spatially separated descendants.

The joint target information remains one bit for exact copies, while the
connected-lesion distance equals the spatial span of the carriers on a path.
"""

from math import ceil
from typing import Tuple


def max_connected_lesion_distance(
    n: int,
    left0: int,
    right0: int,
    t: int,
    speed: int = 1,
) -> int:
    """1-indexed path. Exact maximal span under outward speed bound."""
    if not (1 <= left0 <= right0 <= n):
        raise ValueError("invalid initial carrier interval")
    if t < 0 or speed < 0:
        raise ValueError("t and speed must be nonnegative")

    left_gain = min(left0 - 1, speed * t)
    right_gain = min(n - right0, speed * t)
    initial_span = right0 - left0 + 1
    return initial_span + left_gain + right_gain


def time_to_lesion_tolerance(
    L: int,
    speed: int = 1,
) -> int:
    """Boundary-free lower bound starting from one carrier."""
    if L < 0 or speed <= 0:
        raise ValueError("L >=0 and speed >0 required")
    return ceil(L / (2 * speed))


def demo() -> str:
    n = 101
    c = 51
    lines = [
        "Lineage geometric-robustness amplification",
        "------------------------------------------",
        "A single exact target bit is copied; joint target MI stays 1 bit.",
    ]
    for t in (0, 1, 5, 10, 20, 25):
        d = max_connected_lesion_distance(n, c, c, t, 1)
        lines.append(
            f"t={t:2d}: maximal connected-lesion distance={d}, "
            f"tolerated lesion size={d-1}"
        )
    lines.append(
        f"minimum boundary-free time to tolerate L=20 at speed 1: "
        f"{time_to_lesion_tolerance(20,1)}"
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
