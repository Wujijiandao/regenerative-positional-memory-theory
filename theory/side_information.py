"""Finite binary examples for the side-information regenerative converse."""
from math import log2

def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p*log2(p) - (1-p)*log2(1-p)

def binary_symmetric_conditional_entropy(error):
    """For equiprobable U and S=U xor Bernoulli(error), H(U|S)=h2(error)."""
    if not 0 <= error <= 0.5:
        raise ValueError("error must be in [0,0.5]")
    return h2(error)

def local_zero_error_bits_required(side_error):
    return binary_symmetric_conditional_entropy(side_error)

def demo():
    rows = []
    for q in (0.0, 0.1, 0.25, 0.5):
        rows.append((q, local_zero_error_bits_required(q)))
    return rows

if __name__ == "__main__":
    for q,b in demo():
        print(f"side error={q:.2f}: residual local zero-error entropy={b:.6f} bits")
