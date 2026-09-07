"""Minimal planarian-inspired AP wound-polarity toy.

This is NOT a fitted biological model.

It illustrates a logic distinction:
    generic wound activation != positional address.

Biological inspiration:
- wnt1 is wound induced and promotes posterior Wnt signaling.
- notum is preferentially induced at anterior-facing wounds and inhibits Wnt.
- pre-existing tissue state is needed to polarize the early wound response.

We encode an abstract binary positional cue rather than claiming a complete
molecular model.
"""

from dataclasses import dataclass
from typing import Optional, Dict


@dataclass(frozen=True)
class WoundState:
    wound: int
    positional_cue: Optional[int]  # 1=anterior-biased, 0=posterior-biased, None=erased


def evaluate_wound(state: WoundState) -> Dict[str, Optional[int]]:
    """Evaluate a minimal Boolean wound GRN.

    Rules:
      wnt1 = wound
      notum = wound AND positional_cue
      beta_cat = wnt1 AND NOT notum
      head_selector = notum
      tail_selector = beta_cat

    If positional_cue is erased, notum and downstream fate are unknown.
    """
    wound = int(bool(state.wound))
    if state.positional_cue is None:
        return {
            "wnt1": wound,
            "notum": None,
            "beta_cat": None,
            "head_selector": None,
            "tail_selector": None,
        }

    cue = int(bool(state.positional_cue))
    wnt1 = wound
    notum = wound & cue
    beta_cat = wnt1 & (1 - notum)
    return {
        "wnt1": wnt1,
        "notum": notum,
        "beta_cat": beta_cat,
        "head_selector": notum,
        "tail_selector": beta_cat,
    }


def fragment(asymmetry: bool = True):
    """Return left/right wound states.

    With asymmetry:
        left is anterior-biased; right is posterior-biased.

    Without asymmetry:
        both wounds retain only the same generic wound input and their
        positional cue is erased.
    """
    if asymmetry:
        left = WoundState(1, 1)
        right = WoundState(1, 0)
    else:
        left = WoundState(1, None)
        right = WoundState(1, None)
    return evaluate_wound(left), evaluate_wound(right)


def fate_label(out):
    if out["head_selector"] == 1 and out["tail_selector"] == 0:
        return "HEAD"
    if out["tail_selector"] == 1 and out["head_selector"] == 0:
        return "TAIL"
    return "UNRESOLVED"


def demo():
    a_left, a_right = fragment(asymmetry=True)
    s_left, s_right = fragment(asymmetry=False)

    lines = [
        "Planarian-inspired Boolean GRN toy",
        "----------------------------------",
        f"With positional asymmetry: left={fate_label(a_left)}, right={fate_label(a_right)}",
        f"Generic wnt1 response: left={a_left['wnt1']}, right={a_right['wnt1']}",
        f"With positional cue erased: left={fate_label(s_left)}, right={fate_label(s_right)}",
        f"Generic wound program remains active after cue erasure: "
        f"{s_left['wnt1'] == 1 and s_right['wnt1'] == 1}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(demo())
