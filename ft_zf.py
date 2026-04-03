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

ftz(g, faults=1, return_sets=False)
    Short alias for fault_tolerant_zero_forcing_number.

Z(g, return_sets=False)
    Short alias for zero_forcing_number.

load_all()
    Return a dict mapping every public name in this module to its callable.

Graph formats accepted
----------------------
- NetworkX graph  (has .nodes() and .neighbors(v))
- SageMath graph  (has .vertices() and .neighbors(v))
- dict            {vertex: iterable_of_neighbours}
"""

from itertools import combinations

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

    Examples
    --------
    >>> import networkx as nx
    >>> zero_forcing_number(nx.path_graph(5))
    1
    >>> zero_forcing_number(nx.cycle_graph(4))
    2
    >>> zero_forcing_number(nx.complete_graph(4))
    3
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    if n == 0:
        return 0
    full_mask = (1 << n) - 1
    return _zero_forcing_number_internal(adj_mask, n, full_mask)


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

    Notes
    -----
    **Algorithm overview**

    1. Compute the zero forcing number Z(G) to establish the lower bound
       ftZ(G, k) ≥ Z(G) + k.
    2. Iterate over vertex subsets in increasing size order starting from
       Z(G) + k.
    3. For each candidate set S of size *s*, check whether every
       (s − k)-element subset is a zero forcing set.
    4. Use a shared memoisation cache (mask → zero forcing closure) to avoid
       redundant propagation computations across all subset checks.
    5. Stop as soon as the minimum size is found (with early exit if
       ``return_sets=False``).

    **Complexity**

    In the worst case the algorithm enumerates O(C(n, ftZ(G,k))) candidate
    sets and performs O(C(ftZ, k)) zero forcing checks per candidate, each
    O(n²).  Memoisation of closures substantially reduces repeated work.
    The bitmask representation keeps constants small for graphs up to ~60
    vertices; for larger graphs performance degrades but remains correct.

    Examples
    --------
    >>> import networkx as nx

    Path graph – ftZ equals 2 (the two endpoints):

    >>> fault_tolerant_zero_forcing_number(nx.path_graph(5))
    2

    ``faults=0`` reduces to the standard zero forcing number:

    >>> fault_tolerant_zero_forcing_number(nx.path_graph(5), faults=0)
    1

    Return all minimum sets:

    >>> num, sets = fault_tolerant_zero_forcing_number(
    ...     nx.path_graph(5), faults=1, return_sets=True)
    >>> num
    2
    >>> sorted(sorted(s) for s in sets)
    [[0, 4]]

    Plain dict input:

    >>> g = {0: [1], 1: [0, 2], 2: [1]}   # path P_3
    >>> fault_tolerant_zero_forcing_number(g)
    2
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

def ftz(g, faults: int = 1, return_sets: bool = False):
    """Short alias for :func:`fault_tolerant_zero_forcing_number`."""
    return fault_tolerant_zero_forcing_number(g, faults=faults, return_sets=return_sets)


def Z(g, return_sets: bool = False):
    """Short alias for :func:`zero_forcing_number`.

    When ``return_sets=True`` the behaviour matches
    ``fault_tolerant_zero_forcing_number(g, faults=0, return_sets=True)``.
    """
    return fault_tolerant_zero_forcing_number(g, faults=0, return_sets=return_sets)


def load_all() -> dict:
    """
    Return all public API callables exported by this module.

    Loading is instantaneous (no external files are read) and the returned
    mapping always reflects the live function objects, so the call is safe
    to make at any time and is fully idempotent.

    Returns
    -------
    dict[str, callable]
        A mapping from name to callable for every public function in
        ``ft_zf``:

        - ``"fault_tolerant_zero_forcing_number"`` — compute ftZ(G, k).
        - ``"zero_forcing_number"`` — compute Z(G) (equivalent to faults=0).
        - ``"ftz"`` — short alias for ``fault_tolerant_zero_forcing_number``.
        - ``"Z"`` — short alias for ``zero_forcing_number``.
        - ``"load_all"`` — this function.

    Examples
    --------
    >>> import networkx as nx
    >>> api = load_all()
    >>> api["zero_forcing_number"](nx.path_graph(5))
    1
    >>> api["fault_tolerant_zero_forcing_number"](nx.path_graph(5))
    2
    >>> api["ftz"](nx.path_graph(5), faults=0)
    1
    >>> api2 = load_all()          # idempotent: second call is safe
    >>> api == api2
    True
    """
    return {
        "fault_tolerant_zero_forcing_number": fault_tolerant_zero_forcing_number,
        "zero_forcing_number": zero_forcing_number,
        "ftz": ftz,
        "Z": Z,
        "load_all": load_all,
    }

from itertools import combinations, permutations

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
            U.merge_vertices([u, v])  # merges into one vertex (Sage chooses a representative)

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
