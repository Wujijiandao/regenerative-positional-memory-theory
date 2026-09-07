"""Injury-boundary information bottleneck utilities.

This module implements the finite-state sharp case

    T >= latency - 1 + ceil(K / b)

where:
    K       = target bits required,
    b       = effective target-information channels per update,
    latency = minimum causal steps from boundary to decoder.

It also provides a simple streaming relay construction that attains the bound.

This is an artificial GRN/communication toy, not a quantitative biological model.
"""

from dataclasses import dataclass
from math import ceil
from typing import List, Tuple, Optional


def exact_time_lower_bound(target_bits: int, channels_per_step: int, latency: int) -> int:
    if target_bits < 0:
        raise ValueError("target_bits must be >= 0")
    if channels_per_step <= 0:
        if target_bits == 0:
            return max(0, latency - 1)
        raise ValueError("positive target information cannot cross zero-capacity cut")
    if latency <= 0:
        raise ValueError("latency must be >= 1")
    if target_bits == 0:
        return max(0, latency - 1)
    return latency - 1 + ceil(target_bits / channels_per_step)


def fano_required_bits(num_targets: int, error: float) -> float:
    """Uniform target prior Fano lower bound in bits."""
    import math

    if num_targets < 2:
        return 0.0
    if not (0 <= error < 1):
        raise ValueError("error must be in [0,1)")
    if error == 0:
        return math.log2(num_targets)

    h2 = -error * math.log2(error) - (1 - error) * math.log2(1 - error)
    return max(
        0.0,
        math.log2(num_targets) - h2 - error * math.log2(num_targets - 1)
    )


def approximate_time_lower_bound(
    required_information_bits: float,
    capacity_bits_per_step: float,
    latency: int,
) -> float:
    if required_information_bits < 0:
        raise ValueError("required_information_bits must be >= 0")
    if capacity_bits_per_step <= 0:
        if required_information_bits == 0:
            return float(max(0, latency - 1))
        return float("inf")
    if latency <= 0:
        raise ValueError("latency must be >= 1")
    return latency - 1 + required_information_bits / capacity_bits_per_step


@dataclass
class Packet:
    index: int
    bits: Tuple[int, ...]


def chunk_bits(bits: Tuple[int, ...], width: int) -> List[Packet]:
    if width <= 0:
        raise ValueError("width must be > 0")
    return [
        Packet(i // width, tuple(bits[i:i + width]))
        for i in range(0, len(bits), width)
    ]


def streaming_relay_recovery_time(
    bits: Tuple[int, ...],
    channels_per_step: int,
    latency: int,
) -> Tuple[int, Tuple[int, ...]]:
    """Idealized sharp construction.

    At each update, at most `channels_per_step` target bits cross the boundary.
    A packet then experiences a fixed source-to-decoder latency.

    Packet j (1-indexed injection time) arrives at
        latency + j - 1.

    The decoder concatenates packets in order.
    """
    if any(b not in (0, 1) for b in bits):
        raise ValueError("bits must be binary")
    if latency <= 0:
        raise ValueError("latency must be >= 1")
    if channels_per_step <= 0:
        raise ValueError("channels_per_step must be > 0")

    if not bits:
        return max(0, latency - 1), tuple()

    packets = chunk_bits(bits, channels_per_step)
    arrival_times = [latency + j for j in range(len(packets))]
    final_time = max(arrival_times)
    reconstructed = tuple(bit for packet in packets for bit in packet.bits)
    return final_time, reconstructed


def lesion_core_bound(
    cores: List[Tuple[int, float]],
    capacity_bits_per_step: float,
) -> float:
    """Evaluate max over cores.

    `cores` contains (latency, required_information_bits).
    Returns max latency-1 + info/capacity.
    """
    if not cores:
        return 0.0
    return max(
        approximate_time_lower_bound(info, capacity_bits_per_step, latency)
        for latency, info in cores
    )


def effective_capacity(spatial_capacity: float, grn_capacity: float) -> float:
    if spatial_capacity < 0 or grn_capacity < 0:
        raise ValueError("capacities must be nonnegative")
    return min(spatial_capacity, grn_capacity)


def demo() -> str:
    lines = [
        "Injury-boundary information bottleneck toy",
        "-------------------------------------------",
    ]

    K = 8
    b = 2
    latency = 5
    word = (1, 0, 1, 1, 0, 0, 1, 0)

    lb = exact_time_lower_bound(K, b, latency)
    t, recovered = streaming_relay_recovery_time(word, b, latency)

    lines.append(f"K={K} bits, b={b} bits/step, latency={latency}")
    lines.append(f"lower bound T >= {lb}")
    lines.append(f"streaming construction recovery time = {t}")
    lines.append(f"sharp = {t == lb}")
    lines.append(f"recovered correctly = {recovered == word}")

    c_spatial = 6.0
    c_grn = 2.0
    c_eff = effective_capacity(c_spatial, c_grn)
    lines.append(
        f"spatial capacity={c_spatial}, GRN capacity={c_grn}, "
        f"effective capacity={c_eff}"
    )

    req = fano_required_bits(num_targets=8, error=0.05)
    approx_lb = approximate_time_lower_bound(req, c_eff, latency=3)
    lines.append(
        f"8 equiprobable target classes at 5% error require >= {req:.4f} bits; "
        f"time lower bound at C_eff=2, latency=3 is {approx_lb:.4f}"
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
