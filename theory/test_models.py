import unittest
import math
import networkx as nx

from planarian_ap_grn import fragment, fate_label
from regulatory_cuts import (
    make_fragile_graph,
    make_redundant_graph,
    make_three_path_graph,
    min_target_vertex_cut,
)
from boundary_capacity import (
    exact_time_lower_bound,
    streaming_relay_recovery_time,
    fano_required_bits,
    effective_capacity,
    lesion_core_bound,
)
from connected_lesion import (
    two_marker_code,
    connected_lesion_distance,
    minimum_hamming_distance,
    survives_all_connected_lesions,
    brute_force_distinguishable,
    connected_hull_size,
)
from lineage_information import (
    target_seed_mi,
    target_descendant_joint_mi,
    target_descendant_marginal_sum,
    theoretical_seed_mi,
)
from redundancy_latency import (
    exact_access_formula,
    brute_force_access_radius,
    theorem_optimum,
    brute_force_optimum,
    construct_optimal_pair,
)
from robustness_growth import (
    max_connected_lesion_distance,
    time_to_lesion_tolerance,
)


class TestPlanarianToy(unittest.TestCase):
    def test_positional_asymmetry_selects_opposite_fates(self):
        left, right = fragment(asymmetry=True)
        self.assertEqual(fate_label(left), "HEAD")
        self.assertEqual(fate_label(right), "TAIL")

    def test_generic_wound_program_without_address_is_unresolved(self):
        left, right = fragment(asymmetry=False)
        self.assertEqual(left["wnt1"], 1)
        self.assertEqual(right["wnt1"], 1)
        self.assertEqual(fate_label(left), "UNRESOLVED")
        self.assertEqual(fate_label(right), "UNRESOLVED")


class TestRegulatoryCuts(unittest.TestCase):
    def test_cut_connectivities(self):
        for maker, expected in [
            (make_fragile_graph, 1),
            (make_redundant_graph, 2),
            (make_three_path_graph, 3),
        ]:
            G, S, f = maker()
            k, _ = min_target_vertex_cut(G, S, f)
            self.assertEqual(k, expected)


class TestBoundaryInformationBound(unittest.TestCase):
    def test_sharp_streaming_example(self):
        word = (1, 0, 1, 1, 0, 0, 1, 0)
        bound = exact_time_lower_bound(len(word), 2, 5)
        t, recovered = streaming_relay_recovery_time(word, 2, 5)
        self.assertEqual(bound, 8)
        self.assertEqual(t, bound)
        self.assertEqual(recovered, word)

    def test_effective_capacity(self):
        self.assertEqual(effective_capacity(7.0, 2.5), 2.5)

    def test_fano_zero_error(self):
        self.assertAlmostEqual(fano_required_bits(8, 0.0), 3.0)

    def test_core_optimization(self):
        self.assertAlmostEqual(
            lesion_core_bound([(2, 2.0), (4, 4.0), (6, 0.5)], 2.0),
            5.25,
        )


class TestConnectedLesionTheory(unittest.TestCase):
    def test_spatial_separation(self):
        G = nx.path_graph(9)
        separated = two_marker_code(9, 0, 8)
        clustered = two_marker_code(9, 3, 4)
        self.assertEqual(minimum_hamming_distance(separated), 2)
        self.assertEqual(connected_lesion_distance(G, separated), 9)
        self.assertEqual(connected_lesion_distance(G, clustered), 2)

    def test_theorem_matches_bruteforce(self):
        G = nx.path_graph(7)
        for u in range(7):
            for v in range(u + 1, 7):
                code = two_marker_code(7, u, v)
                d = connected_lesion_distance(G, code)
                for L in range(1, 7):
                    self.assertEqual(
                        d > L,
                        brute_force_distinguishable(G, code, L)
                    )


class TestLineageInformation(unittest.TestCase):
    def test_seed_mi_formula(self):
        for p in (0.0, 0.05, 0.2, 0.4, 0.5):
            self.assertAlmostEqual(
                target_seed_mi(p),
                theoretical_seed_mi(p),
                places=10,
            )

    def test_copying_preserves_joint_information(self):
        p = 0.2
        seed = target_seed_mi(p)
        for n in (1, 2, 4, 8, 32):
            self.assertAlmostEqual(
                target_descendant_joint_mi(p, n),
                seed,
                places=10,
            )
        self.assertGreater(
            target_descendant_marginal_sum(p, 16),
            seed,
        )


class TestRedundancyLatencyTheorem(unittest.TestCase):
    def test_exact_formula_matches_bruteforce(self):
        for n in range(5, 20):
            for L in range(1, min(5, n)):
                for a in range(1, n):
                    for b in range(a + 1, n + 1):
                        if b - a < L:
                            continue
                        self.assertEqual(
                            exact_access_formula(n, a, b, L),
                            brute_force_access_radius(n, a, b, L),
                        )

    def test_closed_form_optimum_exhaustively(self):
        # Exhaustive check across a substantial finite range.
        for L in range(1, 8):
            for n in range(3 * L + 1, 3 * L + 18):
                formula = theorem_optimum(n, L)
                brute, _ = brute_force_optimum(n, L)
                self.assertEqual(formula, brute)


    def test_fixed_separation_frontier_exhaustively_v16(self):
        from redundancy_latency import fixed_separation_optimum, brute_force_fixed_separation
        for L in range(1,7):
            for n in range(3*L+1,3*L+24):
                for s in range(L,n-2*L+2):
                    self.assertEqual(fixed_separation_optimum(n,L,s), brute_force_fixed_separation(n,L,s))

    def test_constructed_pair_attains_theorem(self):
        for L in range(1, 8):
            for n in range(3 * L + 1, 3 * L + 18):
                a, b = construct_optimal_pair(n, L)
                self.assertGreaterEqual(b - a, L)
                self.assertEqual(
                    exact_access_formula(n, a, b, L),
                    theorem_optimum(n, L),
                )


class TestRobustnessGrowth(unittest.TestCase):
    def test_linear_amplification_away_from_boundaries(self):
        n = 101
        c = 51
        for t in range(0, 20):
            self.assertEqual(
                max_connected_lesion_distance(n, c, c, t, 1),
                1 + 2 * t,
            )

    def test_boundary_saturation(self):
        self.assertEqual(
            max_connected_lesion_distance(11, 6, 6, 20, 1),
            11,
        )

    def test_time_to_tolerance(self):
        self.assertEqual(time_to_lesion_tolerance(20, 1), 10)
        self.assertEqual(time_to_lesion_tolerance(21, 2), 6)


# v0.8 multi-source tests
from multisource_single_lesion import (
    theorem_radius as multi_theorem_radius,
    brute_force_optimum as multi_brute_force_optimum,
    construct_optimal as multi_construct_optimal,
    worst_radius_single_cell as multi_worst_radius,
)

class TestMultiSourceSingleLesion(unittest.TestCase):
    def test_exact_formula_small_systems(self):
        for m in range(2, 7):
            for n in range(m, min(17, m + 9)):
                brute, _ = multi_brute_force_optimum(n, m)
                self.assertEqual(brute, multi_theorem_radius(n, m))

    def test_construction_attains_formula(self):
        for m in range(2, 12):
            for n in (m, m + 3, m + 10, m + 25):
                S = multi_construct_optimal(n, m)
                self.assertEqual(len(S), m)
                self.assertEqual(
                    multi_worst_radius(n, S),
                    multi_theorem_radius(n, m),
                )


from side_information import (
    binary_symmetric_conditional_entropy,
    local_zero_error_bits_required,
)

class TestSideInformation(unittest.TestCase):
    def test_perfect_side_information_removes_local_identity_load(self):
        self.assertAlmostEqual(local_zero_error_bits_required(0.0), 0.0)

    def test_target_independent_side_information_recovers_one_bit_load(self):
        self.assertAlmostEqual(local_zero_error_bits_required(0.5), 1.0)

    def test_partial_side_information_reduces_but_does_not_eliminate_load(self):
        q = 0.1
        b = local_zero_error_bits_required(q)
        self.assertGreater(b, 0.0)
        self.assertLess(b, 1.0)


from total_distance_domination_baseline import (
    path_single_cell_optimum,
    path_total_distance_domination_number,
    inverse_consistency,
)

class TestTotalDistanceDominationBaseline(unittest.TestCase):
    def test_known_small_paths(self):
        self.assertEqual(path_total_distance_domination_number(4, 1), 2)
        self.assertEqual(path_total_distance_domination_number(5, 1), 3)
        self.assertEqual(path_total_distance_domination_number(6, 1), 4)

    def test_inverse_grid(self):
        for n in range(2, 80):
            for m in range(2, n + 1):
                self.assertTrue(inverse_consistency(n, m))

    def test_project_example(self):
        self.assertEqual(path_single_cell_optimum(40, 4), 7)


if __name__ == "__main__":
    unittest.main(verbosity=2)
