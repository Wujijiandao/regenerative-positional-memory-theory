"""Small causal-cut utilities for regenerative GRN toy graphs.

The graph nodes may be interpreted as cell-gene-time states.
This code uses NetworkX only for basic directed reachability.

The theoretical point:
If target-correlated sources can influence a fate selector only through a
small node cut, making that cut target-independent destroys target specificity.
"""

from itertools import combinations
from typing import Iterable, Set, Tuple, List
import networkx as nx


def is_cut(G: nx.DiGraph, sources: Iterable[str], sink: str, removed: Set[str]) -> bool:
    H = G.copy()
    H.remove_nodes_from(removed)
    if sink not in H:
        return True
    for s in sources:
        if s in H and nx.has_path(H, s, sink):
            return False
    return True


def min_target_vertex_cut(
    G: nx.DiGraph,
    sources: Iterable[str],
    sink: str,
    perturbable: Iterable[str] = None,
) -> Tuple[int, Set[str]]:
    sources = tuple(sources)
    if perturbable is None:
        perturbable = [n for n in G.nodes if n != sink]
    candidates = [n for n in perturbable if n != sink]

    # If there is already no target-information path.
    if is_cut(G, sources, sink, set()):
        return 0, set()

    for k in range(1, len(candidates) + 1):
        for comb in combinations(candidates, k):
            rem = set(comb)
            if is_cut(G, sources, sink, rem):
                return k, rem
    raise RuntimeError("No finite cut found; check graph/sink protection.")


def make_fragile_graph() -> Tuple[nx.DiGraph, List[str], str]:
    """One target-information route: positional state -> selector -> fate."""
    G = nx.DiGraph()
    edges = [
        ("PCG_state@t0", "notum_selector@t1"),
        ("generic_wound@t0", "notum_selector@t1"),
        ("notum_selector@t1", "AP_fate@t2"),
    ]
    G.add_edges_from(edges)
    return G, ["PCG_state@t0"], "AP_fate@t2"


def make_redundant_graph() -> Tuple[nx.DiGraph, List[str], str]:
    """Two independent target-information routes into a fate integrator.

    Names are abstract: this is not a claim that these exact independent
    pathways are established in planaria.
    """
    G = nx.DiGraph()
    edges = [
        ("positional_module_A@t0", "decoder_A@t1"),
        ("decoder_A@t1", "fate_integrator@t2"),
        ("positional_module_B@t0", "decoder_B@t1"),
        ("decoder_B@t1", "fate_integrator@t2"),
        ("generic_wound@t0", "fate_integrator@t2"),
    ]
    G.add_edges_from(edges)
    return G, ["positional_module_A@t0", "positional_module_B@t0"], "fate_integrator@t2"


def make_three_path_graph() -> Tuple[nx.DiGraph, List[str], str]:
    G = nx.DiGraph()
    for i in range(3):
        s = f"positional_{i}@t0"
        d = f"decoder_{i}@t1"
        G.add_edge(s, d)
        G.add_edge(d, "fate@t2")
    G.add_edge("generic_wound@t0", "fate@t2")
    return G, [f"positional_{i}@t0" for i in range(3)], "fate@t2"


def demo():
    lines = ["Regenerative regulatory connectivity toys", "----------------------------------------"]
    for name, maker in [
        ("fragile", make_fragile_graph),
        ("two-path redundant", make_redundant_graph),
        ("three-path redundant", make_three_path_graph),
    ]:
        G, S, f = maker()
        k, cut = min_target_vertex_cut(G, S, f)
        lines.append(f"{name}: kappa_reg={k}, example_cut={sorted(cut)}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
