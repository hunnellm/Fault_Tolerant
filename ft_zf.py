#!/usr/bin/env python3
"""
ft_zf.py - Fast implementation of the fault tolerant zero forcing number.

A *k*-fault tolerant zero forcing set (k-FTZF set) of a graph G is a vertex
set S such that every subset of S with |S|−k elements is a zero forcing set
of G.  When k = 0 this reduces to the standard zero forcing number Z(G).

The fault tolerant zero forcing number ftZ(G, k) is the minimum size of any
k-FTZF set.

Public API
----------
fault_tolerant_zero_forcing_number(g, faults=1, return_sets=False)
    Compute the k-fault tolerant zero forcing number of a graph.

zero_forcing_number(g)
    Compute the standard zero forcing number (equivalent to faults=0).

propagation_time(g, initial_set)
    Compute the zero forcing propagation time of a given initial set.

direct_fault_propagation_time(g, faults=1, return_sets=False)
    Compute the direct fault propagation time dfpt(G).

indirect_fault_propagation_time_of_set(g, B)
    Compute ifpt(G, B) for a minimum fault tolerant zero forcing set B.

indirect_fault_propagation_time(g, faults=1, return_sets=False)
    Compute the indirect fault propagation time ifpt(G).

ftz(g, faults=1, return_sets=False)
    Short alias for fault_tolerant_zero_forcing_number.

Z(g, return_sets=False)
    Short alias for zero_forcing_number.

dfpt(g, faults=1, return_sets=False)
    Short alias for direct_fault_propagation_time.

ifpt(g, faults=1, return_sets=False)
    Short alias for indirect_fault_propagation_time.

load_all()
    Return a dict mapping every public name in this module to its callable.

Graph formats accepted
----------------------
- NetworkX graph  (has .nodes() and .neighbors(v))
- SageMath graph  (has .vertices() and .neighbors(v))
- dict            {vertex: iterable_of_neighbours}
"""

from itertools import combinations, permutations

# ---------------------------------------------------------------------------
# Graph normalisation
# ---------------------------------------------------------------------------

def _adjacency_lists(g):
    """
    Return ``(vertices, adj_mask, n)`` in a bitmask-ready format.

    Parameters
    ----------
    g : graph or dict
        Any supported graph representation.

    Returns
    -------
    vertices : list
        Sorted list of vertex labels (used to map labels → indices).
    adj_mask : list of int
        ``adj_mask[i]`` is a bitmask of the neighbours of vertex ``i``.
    n : int
        Number of vertices.
    """
    if isinstance(g, dict):
        vertices = sorted(g.keys())
        raw_adj = {v: list(g[v]) for v in vertices}
    elif hasattr(g, 'adjacency'):
        # NetworkX graph
        adj_raw = dict(g.adjacency())
        vertices = sorted(adj_raw.keys())
        raw_adj = {v: list(adj_raw[v].keys()) for v in vertices}
    else:
        # SageMath-compatible (has .vertices() and .neighbors())
        vertices = sorted(g.vertices())
        raw_adj = {v: list(g.neighbors(v)) for v in vertices}

    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}

    adj_mask = [0] * n
    for i, v in enumerate(vertices):
        for u in raw_adj[v]:
            if u in idx:
                j = idx[u]
                adj_mask[i] |= (1 << j)

    return vertices, adj_mask, n


# ---------------------------------------------------------------------------
# Core zero forcing computation (bitmask-based, memoised)
# ---------------------------------------------------------------------------

def _zf_closure(adj_mask, initial_mask, n):
    """
    Compute the zero forcing closure of *initial_mask*.

    Uses the standard colour-change rule: a black vertex v can force its
    unique white neighbour to become black.

    Parameters
    ----------
    adj_mask : list of int
        ``adj_mask[i]`` is a bitmask of the neighbours of vertex ``i``.
    initial_mask : int
        Bitmask of the initially black vertices.
    n : int
        Number of vertices.

    Returns
    -------
    int
        Bitmask of all black vertices after propagation.
    """
    black = initial_mask
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if (black >> v) & 1:          # v is black
                white_nbrs = adj_mask[v] & ~black
                # Check for exactly one white neighbour (non-zero power of two)
                if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                    black |= white_nbrs
                    changed = True
    return black



def _bitmask_from_vertices(vertices, subset):
    """Return the bitmask for ``subset`` using the given ordered ``vertices``."""
    idx = {v: i for i, v in enumerate(vertices)}
    mask = 0
    for v in subset:
        if v not in idx:
            raise ValueError("vertex {!r} is not in the graph".format(v))
        mask |= (1 << idx[v])
    return mask



def _propagation_time_from_mask(adj_mask, initial_mask, n, full_mask):
    """
    Compute propagation time for a zero forcing set given by ``initial_mask``.

    The propagation time is the number of parallel forcing rounds needed to
    color all vertices black, where in each round every currently black vertex
    with a unique white neighbor forces that neighbor.

    Parameters
    ----------
    adj_mask : list of int
        Bitmask adjacency lists.
    initial_mask : int
        Bitmask of initially black vertices.
    n : int
        Number of vertices.
    full_mask : int
        Bitmask with all vertices black.

    Returns
    -------
    int
        The propagation time.

    Raises
    ------
    ValueError
        If ``initial_mask`` is not a zero forcing set.
    """
    black = initial_mask
    steps = 0

    while black != full_mask:
        new_black = 0
        for v in range(n):
            if (black >> v) & 1:
                white_nbrs = adj_mask[v] & ~black
                if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                    new_black |= white_nbrs

        if new_black == 0:
            raise ValueError("initial set is not a zero forcing set")

        black |= new_black
        steps += 1

    return steps


def _zero_forcing_number_internal(adj_mask, n, full_mask):
    """
    Compute Z(G) by testing subsets in increasing size order.

    Returns the size of the smallest zero forcing set found.
    """
    for size in range(0, n + 1):
        for combo in combinations(range(n), size):
            mask = 0
            for v in combo:
                mask |= (1 << v)
            if _zf_closure(adj_mask, mask, n) == full_mask:
                return size
    return n  # pragma: no cover – only reached for empty graph


def _is_ft_zf_set(adj_mask, mask, size, faults, n, full_mask, cache):
    """
    Determine whether *mask* is a k-fault tolerant zero forcing set.

    A set S is k-fault tolerant if every subset of S with |S|−k elements is
    a zero forcing set.  (S itself must therefore also be a ZF set.)

    The *cache* dict (mask → closure) is updated in-place for memoisation.

    Returns
    -------
    bool
    """
    # 1. S itself must be a ZF set
    if mask not in cache:
        cache[mask] = _zf_closure(adj_mask, mask, n)
    if cache[mask] != full_mask:
        return False

    subset_size = size - faults
    if subset_size <= 0:
        # Removing k ≥ |S| vertices leaves the empty set, which can only
        # force a graph with no vertices.
        return False

    # 2. Every (size − faults)-element subset must also be a ZF set
    vertices_in_s = [v for v in range(n) if (mask >> v) & 1]
    for sub in combinations(vertices_in_s, subset_size):
        sub_mask = 0
        for v in sub:
            sub_mask |= (1 << v)
        if sub_mask not in cache:
            cache[sub_mask] = _zf_closure(adj_mask, sub_mask, n)
        if cache[sub_mask] != full_mask:
            return False

    return True


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def zero_forcing_number(g):
    """
    Compute the standard zero forcing number Z(G).

    Equivalent to ``fault_tolerant_zero_forcing_number(g, faults=0)``.

    Parameters
    ----------
    g : graph
        A graph in any supported format (NetworkX, SageMath, or adjacency
        dict ``{v: [neighbours]}``.

    Returns
    -------
    int
        The zero forcing number Z(G).
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    if n == 0:
        return 0
    full_mask = (1 << n) - 1
    return _zero_forcing_number_internal(adj_mask, n, full_mask)



def propagation_time(g, initial_set):
    """
    Compute the zero forcing propagation time of ``initial_set`` on ``g``.

    Parameters
    ----------
    g : graph
        A graph in any format supported by this module.
    initial_set : iterable
        Vertices initially colored black.

    Returns
    -------
    int
        Number of parallel forcing rounds needed to color all vertices.

    Raises
    ------
    ValueError
        If ``initial_set`` is not a zero forcing set.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    full_mask = (1 << n) - 1
    mask = _bitmask_from_vertices(vertices, initial_set)
    return _propagation_time_from_mask(adj_mask, mask, n, full_mask)



def fault_tolerant_zero_forcing_number(g, faults=1, return_sets=False):
    """
    Compute the fault tolerant zero forcing number ftZ(G, faults).

    The *k*-fault tolerant zero forcing number ftZ(G, k) is the size of the
    smallest vertex set S such that every subset of S with |S|−k elements is
    a zero forcing set of G.  Setting ``faults=0`` gives the standard zero
    forcing number Z(G).

    Parameters
    ----------
    g : graph
        A graph in any supported format:

        - **NetworkX** graph (``networkx.Graph`` or ``DiGraph``).
        - **SageMath** graph (objects exposing ``.vertices()`` /
          ``.neighbors()``).
        - **dict** ``{vertex: iterable_of_neighbours}``.

    faults : int, optional
        Number of faults to tolerate (default: **1**).  Must be ≥ 0.

    return_sets : bool, optional
        If ``True``, return a tuple ``(ftZ, sets)`` where *sets* is a sorted
        list of :class:`frozenset` objects, one per minimum k-fault tolerant
        zero forcing set.  Default: **False**.

    Returns
    -------
    int
        The fault tolerant zero forcing number (when ``return_sets=False``).
    tuple ``(int, list of frozenset)``
        A pair ``(ftZ, sets)`` (when ``return_sets=True``).

    Raises
    ------
    ValueError
        If *faults* < 0.
    """
    if faults < 0:
        raise ValueError("faults must be >= 0")

    vertices, adj_mask, n = _adjacency_lists(g)

    if n == 0:
        if return_sets:
            return 0, [frozenset()]
        return 0

    full_mask = (1 << n) - 1

    # ----------------------------------------------------------------
    # faults == 0: reduce to standard zero forcing
    # ----------------------------------------------------------------
    if faults == 0:
        z = _zero_forcing_number_internal(adj_mask, n, full_mask)
        if not return_sets:
            return z
        cache = {}
        zf_sets = []
        for combo in combinations(range(n), z):
            mask = 0
            for v in combo:
                mask |= (1 << v)
            if mask not in cache:
                cache[mask] = _zf_closure(adj_mask, mask, n)
            if cache[mask] == full_mask:
                zf_sets.append(frozenset(vertices[v] for v in combo))
        return z, sorted(zf_sets, key=lambda s: sorted(s))

    # ----------------------------------------------------------------
    # General case: faults >= 1
    # ----------------------------------------------------------------

    # Lower bound: ftZ(G, k) >= Z(G) + k
    z = _zero_forcing_number_internal(adj_mask, n, full_mask)
    lower = z + faults

    cache = {}   # shared memoisation: mask -> ZF closure mask

    ftz = -1
    ftz_sets = []

    for size in range(lower, n + 1):
        if ftz != -1:
            break  # minimum already found; all sets of this size collected

        found_at_this_size = False
        for combo in combinations(range(n), size):
            mask = 0
            for v in combo:
                mask |= (1 << v)

            if _is_ft_zf_set(adj_mask, mask, size, faults, n, full_mask, cache):
                ftz = size
                found_at_this_size = True
                if return_sets:
                    ftz_sets.append(frozenset(vertices[v] for v in combo))
                else:
                    break  # one set suffices when we only need the number

        if found_at_this_size and not return_sets:
            break

    if return_sets:
        return ftz, sorted(ftz_sets, key=lambda s: sorted(s))
    return ftz



def minimum_fault_tolerant_zero_forcing_sets(g, faults=1):
    """
    Return all minimum fault tolerant zero forcing sets of ``g``.

    Parameters
    ----------
    g : graph
        Graph in any format supported by this module.
    faults : int, optional
        Number of tolerated faults. Default is 1.

    Returns
    -------
    tuple
        ``(ftz, sets)`` where ``ftz`` is the minimum size and ``sets`` is the
        list of all minimum fault tolerant zero forcing sets.
    """
    return fault_tolerant_zero_forcing_number(g, faults=faults, return_sets=True)



def direct_fault_propagation_time(g, faults=1, return_sets=False):
    """
    Compute the direct fault propagation time ``dfpt(G)``.

    This is the minimum propagation time among all minimum ``faults``-fault
    tolerant zero forcing sets of ``g``.

    Parameters
    ----------
    g : graph
        Graph in any format supported by this module.
    faults : int, optional
        Number of tolerated faults. Default is 1.
    return_sets : bool, optional
        If True, also return the minimum FTZF sets achieving ``dfpt``.

    Returns
    -------
    int
        The direct fault propagation time.
    tuple
        If ``return_sets=True``, returns ``(dfpt, sets)``.
    """
    _, ft_sets = minimum_fault_tolerant_zero_forcing_sets(g, faults=faults)

    best_time = None
    best_sets = []

    for B in ft_sets:
        t = propagation_time(g, B)
        if best_time is None or t < best_time:
            best_time = t
            best_sets = [B]
        elif t == best_time:
            best_sets.append(B)

    if return_sets:
        return best_time, sorted(best_sets, key=lambda s: sorted(s))
    return best_time



def indirect_fault_propagation_time_of_set(g, B):
    """
    Compute ``ifpt(G, B)`` for a minimum fault tolerant zero forcing set ``B``.

    For each vertex ``v`` in ``B``, compute the propagation time of ``B - {v}``.
    The result is the maximum of those propagation times.

    Parameters
    ----------
    g : graph
        Graph in any format supported by this module.
    B : iterable
        A minimum fault tolerant zero forcing set.

    Returns
    -------
    int
        ``ifpt(G, B)``.
    """
    B = frozenset(B)
    worst_time = 0
    for v in B:
        t = propagation_time(g, B - {v})
        if t > worst_time:
            worst_time = t
    return worst_time



def indirect_fault_propagation_time(g, faults=1, return_sets=False):
    """
    Compute the indirect fault propagation time ``ifpt(G)``.

    This is the minimum value of ``ifpt(G, B)`` among all minimum
    ``faults``-fault tolerant zero forcing sets ``B`` of ``g``.

    Parameters
    ----------
    g : graph
        Graph in any format supported by this module.
    faults : int, optional
        Number of tolerated faults. Default is 1.
    return_sets : bool, optional
        If True, also return the minimum FTZF sets achieving ``ifpt``.

    Returns
    -------
    int
        The indirect fault propagation time.
    tuple
        If ``return_sets=True``, returns ``(ifpt, sets)``.
    """
    _, ft_sets = minimum_fault_tolerant_zero_forcing_sets(g, faults=faults)

    best_time = None
    best_sets = []

    for B in ft_sets:
        t = indirect_fault_propagation_time_of_set(g, B)
        if best_time is None or t < best_time:
            best_time = t
            best_sets = [B]
        elif t == best_time:
            best_sets.append(B)

    if return_sets:
        return best_time, sorted(best_sets, key=lambda s: sorted(s))
    return best_time



def ftz(g, faults=1, return_sets=False):
    """Short alias for :func:`fault_tolerant_zero_forcing_number`."""
    return fault_tolerant_zero_forcing_number(g, faults=faults, return_sets=return_sets)



def Z(g, return_sets=False):
    """Short alias for :func:`zero_forcing_number`.

    When ``return_sets=True`` the behaviour matches
    ``fault_tolerant_zero_forcing_number(g, faults=0, return_sets=True)``.
    """
    return fault_tolerant_zero_forcing_number(g, faults=0, return_sets=return_sets)



def dfpt(g, faults=1, return_sets=False):
    """Short alias for :func:`direct_fault_propagation_time`."""
    return direct_fault_propagation_time(g, faults=faults, return_sets=return_sets)



def ifpt(g, faults=1, return_sets=False):
    """Short alias for :func:`indirect_fault_propagation_time`."""
    return indirect_fault_propagation_time(g, faults=faults, return_sets=return_sets)



def load_all():
    """
    Return all public API callables exported by this module.

    Loading is instantaneous (no external files are read) and the returned
    mapping always reflects the live function objects, so the call is safe
    to make at any time and is fully idempotent.

    Returns
    -------
    dict
        A mapping from name to callable for every public function in
        ``ft_zf``.
    """
    return {
        "fault_tolerant_zero_forcing_number": fault_tolerant_zero_forcing_number,
        "zero_forcing_number": zero_forcing_number,
        "propagation_time": propagation_time,
        "minimum_fault_tolerant_zero_forcing_sets": minimum_fault_tolerant_zero_forcing_sets,
        "direct_fault_propagation_time": direct_fault_propagation_time,
        "indirect_fault_propagation_time_of_set": indirect_fault_propagation_time_of_set,
        "indirect_fault_propagation_time": indirect_fault_propagation_time,
        "ftz": ftz,
        "Z": Z,
        "dfpt": dfpt,
        "ifpt": ifpt,
        "load_all": load_all,
    }


def identify_graphs(g1, g2, identify):
    """
    Return a list of all (non-isomorphic) graphs obtained by identifying
    'identify' vertices of g1 with 'identify' vertices of g2 pairwise.

    INPUT:
      - g1, g2: Sage Graph objects
      - identify: nonnegative integer k

    OUTPUT:
      - list of Sage Graph objects (deduplicated up to isomorphism)

    Notes:
      - Vertices are treated as distinct across g1 and g2 even if they have
        equal labels; the function relabels internally to avoid collisions.
      - The k identifications are done in all possible bijections between
        chosen k-subsets.
    """
    k = int(identify)
    if k < 0:
        raise ValueError("identify must be >= 0")
    if k > g1.num_verts() or k > g2.num_verts():
        return []

    # Relabel so vertex sets are disjoint and stable
    G1 = g1.copy(immutable=False)
    G2 = g2.copy(immutable=False)
    G1.relabel({v: ("g1", v) for v in G1.vertices()}, inplace=True)
    G2.relabel({v: ("g2", v) for v in G2.vertices()}, inplace=True)

    V1 = list(G1.vertices())
    V2 = list(G2.vertices())

    results = []
    canon_seen = set()

    # helper: union graph, then contract vertices according to pairs
    def glue_once(pairs):
        # disjoint union: keep all edges, no edges between components
        U = Graph()
        U.add_vertices(V1)
        U.add_edges(G1.edges(labels=False))
        U.add_vertices(V2)
        U.add_edges(G2.edges(labels=False))

        # Contract each pair (u in g1, v in g2) into u, removing v
        # Important: after contractions, later vertex names still exist,
        # but only if not previously removed.
        for (u, v) in pairs:
            if u == v:
                continue
            if (u not in U) or (v not in U):
                continue
            U.merge_vertices([u, v])  # Sage chooses a representative

        return U

    for A in combinations(V1, k):
        for B in combinations(V2, k):
            for permB in permutations(B):
                pairs = list(zip(A, permB))
                H = glue_once(pairs)

                # Deduplicate up to isomorphism via canonical label
                C = H.canonical_label()
                key = (C.graph6_string() if hasattr(C, "graph6_string") else C.certificate())
                if key not in canon_seen:
                    canon_seen.add(key)
                    results.append(H)

    return results
