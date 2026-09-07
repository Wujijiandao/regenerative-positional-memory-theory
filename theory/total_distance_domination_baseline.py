"""Single-cell access radius as inverse total distance-r domination."""
from math import ceil

def path_single_cell_optimum(n, m):
    if not (n >= m >= 2):
        raise ValueError("require n >= m >= 2")
    c = m // 2
    return max(1, ceil((n - c) / (m + c)))

def path_total_distance_domination_number(n, r):
    if n < 2 or r < 1:
        raise ValueError("require n >= 2 and r >= 1")
    q, s = divmod(n, 3*r + 1)
    if s == 0:
        g = 2*q
    elif s <= r:
        g = 2*q + 1
    else:
        g = 2*q + 2
    return max(2, g)

def inverse_consistency(n, m):
    r = path_single_cell_optimum(n, m)
    return (
        path_total_distance_domination_number(n, r) <= m
        and (r == 1 or path_total_distance_domination_number(n, r-1) > m)
    )
