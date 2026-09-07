"""Connected-lesion distance for spatial regulatory-memory codebooks.

This is a small exact tool for finite tissue graphs.

Definitions
-----------
difference support:
    vertices where two target states differ.

connected hull:
    smallest connected vertex set containing a difference support.

connected-lesion distance:
    minimum connected-hull size over target-state pairs.

For arbitrary erasures, Hamming distance controls zero-error erasure
distinguishability. For connected lesions, geometry can make the effective
distance much larger.
"""

from itertools import combinations
from typing import Dict, Hashable, Iterable, Sequence, Set, Tuple
import networkx as nx


State = Tuple[Hashable, ...]


def difference_support(a: State, b: State) -> Set[int]:
    if len(a) != len(b):
        raise ValueError("states must have equal length")
    return {i for i, (x, y) in enumerate(zip(a, b)) if x != y}


def hamming_distance(a: State, b: State) -> int:
    return len(difference_support(a, b))


def connected_hull_size(G: nx.Graph, support: Iterable[int]) -> int:
    """Exact minimum number of vertices in a connected set containing support.

    Exponential brute force; intended for small theorem toys.
    """
    S = set(support)
    if not S:
        return 0
    if not S.issubset(G.nodes):
        raise ValueError("support contains vertices not in graph")

    nodes = list(G.nodes)
    extras = [v for v in nodes if v not in S]

    for k in range(len(S), len(nodes) + 1):
        need = k - len(S)
        for add in combinations(extras, need):
            D = S.union(add)
            if nx.is_connected(G.subgraph(D)):
                return k
    raise RuntimeError("No connected hull found; graph may be disconnected")


def path_hull_size(support: Iterable[int]) -> int:
    S = list(support)
    if not S:
        return 0
    return max(S) - min(S) + 1


def connected_lesion_distance(
    G: nx.Graph,
    codebook: Dict[Hashable, State],
) -> int:
    labels = list(codebook)
    if len(labels) < 2:
        raise ValueError("need at least two target states")

    best = None
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a = codebook[labels[i]]
            b = codebook[labels[j]]
            S = difference_support(a, b)
            if not S:
                return 0
            h = connected_hull_size(G, S)
            best = h if best is None else min(best, h)
    return best


def minimum_hamming_distance(codebook: Dict[Hashable, State]) -> int:
    labels = list(codebook)
    if len(labels) < 2:
        raise ValueError("need at least two target states")
    return min(
        hamming_distance(codebook[labels[i]], codebook[labels[j]])
        for i in range(len(labels))
        for j in range(i + 1, len(labels))
    )


def survives_all_connected_lesions(
    G: nx.Graph,
    codebook: Dict[Hashable, State],
    max_lesion_size: int,
) -> bool:
    return connected_lesion_distance(G, codebook) > max_lesion_size


def enumerate_connected_sets(G: nx.Graph, max_size: int):
    nodes = list(G.nodes)
    for k in range(1, max_size + 1):
        for D in combinations(nodes, k):
            if nx.is_connected(G.subgraph(D)):
                yield set(D)


def erased_observation(state: State, D: Set[int]):
    return tuple(None if i in D else x for i, x in enumerate(state))


def brute_force_distinguishable(
    G: nx.Graph,
    codebook: Dict[Hashable, State],
    max_lesion_size: int,
) -> bool:
    labels = list(codebook)
    for D in enumerate_connected_sets(G, max_lesion_size):
        seen = set()
        for lab in labels:
            obs = erased_observation(codebook[lab], D)
            if obs in seen:
                return False
            seen.add(obs)
    return True


def two_marker_code(n: int, u: int, v: int):
    """Two binary targets differing only at u and v."""
    if not (0 <= u < n and 0 <= v < n and u != v):
        raise ValueError("invalid marker locations")
    a = [0] * n
    b = [0] * n
    b[u] = 1
    b[v] = 1
    return {"A": tuple(a), "B": tuple(b)}


def demo() -> str:
    n = 9
    G = nx.path_graph(n)

    separated = two_marker_code(n, 0, 8)
    clustered = two_marker_code(n, 3, 4)

    lines = [
        "Connected-lesion positional-memory toy",
        "--------------------------------------",
        f"path length n={n}",
        "",
        "Separated two-marker code:",
        f"  Hamming distance = {minimum_hamming_distance(separated)}",
        f"  connected-lesion distance = {connected_lesion_distance(G, separated)}",
        f"  survives every connected lesion <= 8 cells = "
        f"{survives_all_connected_lesions(G, separated, 8)}",
        "",
        "Clustered two-marker code:",
        f"  Hamming distance = {minimum_hamming_distance(clustered)}",
        f"  connected-lesion distance = {connected_lesion_distance(G, clustered)}",
        f"  survives every connected lesion <= 2 cells = "
        f"{survives_all_connected_lesions(G, clustered, 2)}",
        f"  survives every connected lesion <= 1 cell = "
        f"{survives_all_connected_lesions(G, clustered, 1)}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
