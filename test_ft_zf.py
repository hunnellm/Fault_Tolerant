#!/usr/bin/env python3
"""
test_ft_zf.py - Tests for the fault tolerant zero forcing number implementation.

Run with:  pytest test_ft_zf.py -v
"""

import pytest
import networkx as nx

from ft_zf import fault_tolerant_zero_forcing_number, zero_forcing_number, load_all


# ---------------------------------------------------------------------------
# Helper: simple graph builder using dict representation
# ---------------------------------------------------------------------------

def path_dict(n):
    """Return adjacency dict for path P_n (vertices 0..n-1)."""
    if n == 0:
        return {}
    if n == 1:
        return {0: []}
    adj = {}
    for v in range(n):
        nbrs = []
        if v > 0:
            nbrs.append(v - 1)
        if v < n - 1:
            nbrs.append(v + 1)
        adj[v] = nbrs
    return adj


def cycle_dict(n):
    """Return adjacency dict for cycle C_n (vertices 0..n-1)."""
    adj = {}
    for v in range(n):
        adj[v] = [(v - 1) % n, (v + 1) % n]
    return adj


def complete_dict(n):
    """Return adjacency dict for complete graph K_n (vertices 0..n-1)."""
    adj = {}
    for v in range(n):
        adj[v] = [u for u in range(n) if u != v]
    return adj


# ---------------------------------------------------------------------------
# zero_forcing_number  (faults=0 baseline)
# ---------------------------------------------------------------------------

class TestZeroForcingNumber:
    """Tests for the zero_forcing_number helper function."""

    def test_path_p2(self):
        assert zero_forcing_number(nx.path_graph(2)) == 1

    def test_path_p5(self):
        # Z(P_5) = 1  (an endpoint forces its way across)
        assert zero_forcing_number(nx.path_graph(5)) == 1

    def test_cycle_c4(self):
        # Z(C_4) = 2
        assert zero_forcing_number(nx.cycle_graph(4)) == 2

    def test_cycle_c5(self):
        # Z(C_5) = 2
        assert zero_forcing_number(nx.cycle_graph(5)) == 2

    def test_complete_k4(self):
        # Z(K_4) = 3  (need n-1 vertices for a complete graph)
        assert zero_forcing_number(nx.complete_graph(4)) == 3

    def test_complete_k1(self):
        # Z(K_1) = 1
        assert zero_forcing_number(nx.complete_graph(1)) == 1

    def test_dict_input(self):
        # Same as P_5 but given as a dict
        g = path_dict(5)
        assert zero_forcing_number(g) == 1

    def test_empty_graph_zero_vertices(self):
        assert zero_forcing_number({}) == 0

    def test_single_vertex_no_edges(self):
        # A single isolated vertex: {0} forces itself
        g = {0: []}
        assert zero_forcing_number(g) == 1


# ---------------------------------------------------------------------------
# fault_tolerant_zero_forcing_number – basic correctness
# ---------------------------------------------------------------------------

class TestFaultTolerantBasic:
    """Correctness tests for fault_tolerant_zero_forcing_number."""

    # ------------------------------------------------------------------
    # faults=0 must equal Z(G)
    # ------------------------------------------------------------------

    def test_faults0_path_p5(self):
        g = nx.path_graph(5)
        assert fault_tolerant_zero_forcing_number(g, faults=0) == zero_forcing_number(g)

    def test_faults0_cycle_c4(self):
        g = nx.cycle_graph(4)
        assert fault_tolerant_zero_forcing_number(g, faults=0) == zero_forcing_number(g)

    def test_faults0_complete_k4(self):
        g = nx.complete_graph(4)
        assert fault_tolerant_zero_forcing_number(g, faults=0) == zero_forcing_number(g)

    # ------------------------------------------------------------------
    # Default faults=1
    # ------------------------------------------------------------------

    def test_default_faults_is_one(self):
        """Calling without faults kwarg uses faults=1."""
        g = nx.path_graph(5)
        assert fault_tolerant_zero_forcing_number(g) == \
               fault_tolerant_zero_forcing_number(g, faults=1)

    def test_path_p5_faults1(self):
        # ftZ(P_5, 1) = 2  (the two endpoints {0, 4})
        assert fault_tolerant_zero_forcing_number(nx.path_graph(5), faults=1) == 2

    def test_path_p3_faults1(self):
        # ftZ(P_3, 1) = 2  ({0, 2}: each endpoint is a ZF set)
        assert fault_tolerant_zero_forcing_number(nx.path_graph(3), faults=1) == 2

    def test_path_p2_faults1(self):
        # ftZ(P_2, 1) = 2  (both singletons are ZF sets)
        assert fault_tolerant_zero_forcing_number(nx.path_graph(2), faults=1) == 2

    def test_complete_k4_faults1(self):
        # ftZ(K_4, 1) = 4  (all vertices needed)
        assert fault_tolerant_zero_forcing_number(nx.complete_graph(4), faults=1) == 4

    def test_complete_k3_faults1(self):
        # Z(K_3) = 2; ftZ(K_3, 1) = 3
        assert fault_tolerant_zero_forcing_number(nx.complete_graph(3), faults=1) == 3

    def test_cycle_c4_faults1(self):
        # Z(C_4) = 2; ftZ(C_4, 1) = 4
        # (C_4 has no triangle, so no 3-element set has all adjacent pairs)
        assert fault_tolerant_zero_forcing_number(nx.cycle_graph(4), faults=1) == 4

    def test_cycle_c5_faults1(self):
        # Z(C_5) = 2; ftZ(C_5, 1) = 4
        assert fault_tolerant_zero_forcing_number(nx.cycle_graph(5), faults=1) == 4

    def test_path_p5_faults2(self):
        # ftZ(P_5, 2) = 4.
        # Only {0} and {4} are ZF singletons, so no 3-element set can have
        # all its singletons be ZF sets.  The minimum 2-FT ZF sets of size 4
        # are {0,1,2,4} and {0,2,3,4}.
        result = fault_tolerant_zero_forcing_number(nx.path_graph(5), faults=2)
        assert result == 4

    # ------------------------------------------------------------------
    # Dict input
    # ------------------------------------------------------------------

    def test_dict_path_p5_faults1(self):
        g = path_dict(5)
        assert fault_tolerant_zero_forcing_number(g, faults=1) == 2

    def test_dict_complete_k4_faults1(self):
        g = complete_dict(4)
        assert fault_tolerant_zero_forcing_number(g, faults=1) == 4

    # ------------------------------------------------------------------
    # Edge cases
    # ------------------------------------------------------------------

    def test_empty_graph(self):
        assert fault_tolerant_zero_forcing_number({}, faults=1) == 0

    def test_negative_faults_raises(self):
        with pytest.raises(ValueError):
            fault_tolerant_zero_forcing_number(nx.path_graph(3), faults=-1)


# ---------------------------------------------------------------------------
# return_sets=False vs return_sets=True
# ---------------------------------------------------------------------------

class TestReturnSets:
    """Tests for the return_sets parameter."""

    def test_return_sets_false_returns_int(self):
        result = fault_tolerant_zero_forcing_number(
            nx.path_graph(5), faults=1, return_sets=False)
        assert isinstance(result, int)

    def test_return_sets_true_returns_tuple(self):
        result = fault_tolerant_zero_forcing_number(
            nx.path_graph(5), faults=1, return_sets=True)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_return_sets_true_number_matches(self):
        g = nx.path_graph(5)
        num_only = fault_tolerant_zero_forcing_number(g, faults=1)
        num_with_sets, sets = fault_tolerant_zero_forcing_number(
            g, faults=1, return_sets=True)
        assert num_only == num_with_sets

    def test_return_sets_true_sets_are_frozensets(self):
        _, sets = fault_tolerant_zero_forcing_number(
            nx.path_graph(5), faults=1, return_sets=True)
        assert all(isinstance(s, frozenset) for s in sets)

    def test_return_sets_correct_path_p5_faults1(self):
        num, sets = fault_tolerant_zero_forcing_number(
            nx.path_graph(5), faults=1, return_sets=True)
        assert num == 2
        # The only minimum 1-FT ZF set for P_5 is {0, 4}
        assert frozenset({0, 4}) in sets

    def test_return_sets_all_valid_ft_sets(self):
        """Every returned set must actually be a k-FT ZF set."""
        g = nx.path_graph(5)
        num, sets = fault_tolerant_zero_forcing_number(
            g, faults=1, return_sets=True)
        vertices, adj_mask, n = _unpack(g)
        full_mask = (1 << n) - 1
        idx = {v: i for i, v in enumerate(vertices)}
        for s in sets:
            assert len(s) == num
            mask = sum(1 << idx[v] for v in s)
            # Check s is ZF
            from ft_zf import _zf_closure
            assert _zf_closure(adj_mask, mask, n) == full_mask
            # Check every (|s|-1) subset is ZF
            for sub in _subsets_minus_one(s):
                sub_mask = sum(1 << idx[v] for v in sub)
                assert _zf_closure(adj_mask, sub_mask, n) == full_mask

    def test_return_sets_no_duplicates(self):
        _, sets = fault_tolerant_zero_forcing_number(
            nx.cycle_graph(5), faults=1, return_sets=True)
        assert len(sets) == len(set(sets))

    def test_return_sets_faults0_matches_zero_forcing_sets(self):
        """With faults=0, return_sets=True must return all minimum ZF sets."""
        g = nx.path_graph(3)
        num, sets = fault_tolerant_zero_forcing_number(g, faults=0, return_sets=True)
        assert num == 1  # Z(P_3) = 1
        # P_3 has two ZF sets of size 1: {0} and {2}
        assert frozenset({0}) in sets
        assert frozenset({2}) in sets

    def test_return_sets_complete_k4_faults1(self):
        num, sets = fault_tolerant_zero_forcing_number(
            nx.complete_graph(4), faults=1, return_sets=True)
        assert num == 4
        # Only the full vertex set is a 1-FT ZF set for K_4
        assert frozenset({0, 1, 2, 3}) in sets

    def test_return_sets_sorted_deterministic(self):
        """Two identical calls must return the same ordered list of sets."""
        g = nx.cycle_graph(5)
        _, s1 = fault_tolerant_zero_forcing_number(g, faults=1, return_sets=True)
        _, s2 = fault_tolerant_zero_forcing_number(g, faults=1, return_sets=True)
        assert s1 == s2


# ---------------------------------------------------------------------------
# Monotonicity properties
# ---------------------------------------------------------------------------

class TestMonotonicity:
    """Structural sanity checks."""

    def test_ftZ_geq_Z_plus_faults(self):
        """ftZ(G, k) >= Z(G) + k when a k-FT ZF set exists.

        -1 is returned when no k-FT ZF set exists for the given n and k
        (e.g. P_2 with faults=2, where faults >= n leaves only the empty
        subset which cannot force a non-empty graph).
        """
        for n in range(2, 6):
            g = nx.path_graph(n)
            z = zero_forcing_number(g)
            for k in range(1, 3):
                ftz = fault_tolerant_zero_forcing_number(g, faults=k)
                if ftz != -1:   # -1 means no k-FT ZF set exists
                    assert ftz >= z + k

    def test_ftZ_nondecreasing_in_faults(self):
        """ftZ(G, k) <= ftZ(G, k+1)."""
        g = nx.path_graph(6)
        prev = fault_tolerant_zero_forcing_number(g, faults=0)
        for k in range(1, 4):
            curr = fault_tolerant_zero_forcing_number(g, faults=k)
            assert curr >= prev
            prev = curr

    def test_ftZ_leq_n(self):
        """ftZ(G, k) <= n (full set is always k-FT)."""
        for n in range(2, 6):
            g = nx.path_graph(n)
            for k in range(1, 3):
                ftz = fault_tolerant_zero_forcing_number(g, faults=k)
                assert ftz <= n


# ---------------------------------------------------------------------------
# Internal helpers used by tests
# ---------------------------------------------------------------------------

def _unpack(g):
    """Return (vertices, adj_mask, n) for a NetworkX graph."""
    from ft_zf import _adjacency_lists
    return _adjacency_lists(g)


def _subsets_minus_one(s):
    """Yield all subsets of frozenset s with one element removed."""
    lst = sorted(s)
    for i in range(len(lst)):
        yield frozenset(lst[:i] + lst[i + 1:])


# ---------------------------------------------------------------------------
# load_all
# ---------------------------------------------------------------------------

class TestLoadAll:
    """Tests for the load_all() convenience function."""

    EXPECTED_KEYS = {
        "fault_tolerant_zero_forcing_number",
        "zero_forcing_number",
        "ftz",
        "Z",
        "load_all",
    }

    def test_returns_dict(self):
        result = load_all()
        assert isinstance(result, dict)

    def test_contains_expected_keys(self):
        result = load_all()
        assert self.EXPECTED_KEYS == set(result.keys())

    def test_all_values_are_callable(self):
        for name, fn in load_all().items():
            assert callable(fn), f"'{name}' is not callable"

    def test_idempotent(self):
        """Calling load_all() twice returns equivalent mappings."""
        assert load_all() == load_all()

    def test_zero_forcing_number_callable(self):
        api = load_all()
        assert api["zero_forcing_number"](nx.path_graph(5)) == 1

    def test_fault_tolerant_zero_forcing_number_callable(self):
        api = load_all()
        assert api["fault_tolerant_zero_forcing_number"](nx.path_graph(5)) == 2

    def test_ftz_alias_callable(self):
        api = load_all()
        assert api["ftz"](nx.path_graph(5), faults=1) == 2

    def test_ftz_alias_respects_faults_param(self):
        api = load_all()
        # faults=0 should equal the standard zero forcing number
        assert api["ftz"](nx.path_graph(5), faults=0) == 1

    def test_Z_alias_callable(self):
        api = load_all()
        assert api["Z"](nx.path_graph(5)) == 1

    def test_Z_alias_return_sets(self):
        api = load_all()
        num, sets = api["Z"](nx.path_graph(3), return_sets=True)
        assert num == 1
        assert frozenset({0}) in sets
        assert frozenset({2}) in sets

    def test_functions_are_live_references(self):
        """Values in the dict must be the same objects as the module-level names."""
        import ft_zf
        api = load_all()
        assert api["fault_tolerant_zero_forcing_number"] is ft_zf.fault_tolerant_zero_forcing_number
        assert api["zero_forcing_number"] is ft_zf.zero_forcing_number
        assert api["load_all"] is ft_zf.load_all
