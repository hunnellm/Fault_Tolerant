#!/usr/bin/env python3
"""
loop_zf.py - Looped graph zero forcing number for a fixed loop configuration.

A loop configuration of a graph G is obtained by placing loops on a chosen
subset of vertices.  Under the looped color-change rule:

    If any vertex has exactly one white neighbor (possibly itself),
    then that neighbor can be colored blue.

The loop graph zero forcing number for a fixed loop configuration is the
minimum size of an initial blue set that eventually forces all vertices blue.

Public API
----------
looped_zero_forcing_number(g, looped_vertices=None, return_sets=False)
    Compute the looped zero forcing number for the specified loop set.

looped_zero_forcing_closure(g, initial_set, looped_vertices=None)
    Compute the final blue closure of an initial set under the looped rule.

is_looped_zero_forcing_set(g, initial_set, looped_vertices=None)
    Test whether an initial set is looped zero forcing.

lzf(g, looped_vertices=None, return_sets=False)
    Short alias for looped_zero_forcing_number.

load_all()
    Return a dict mapping every public name in this module to its callable.

Notes
-----
- This implementation is designed to work in SageMath and mirrors the style
  used in ft_zf.py (bitmask propagation + subset search in increasing size).
- The input graph g is treated as the underlying simple graph.  Loops are
  specified separately via `looped_vertices`.
"""

from itertools import combinations


# ---------------------------------------------------------------------------
# Graph normalization
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
    elif hasattr(g, "adjacency"):
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


def _bitmask_from_vertices(vertices, subset):
    """Return the bitmask for ``subset`` using the given ordered ``vertices``."""
    idx = {v: i for i, v in enumerate(vertices)}
    mask = 0
    for v in subset:
        if v not in idx:
            raise ValueError("vertex {!r} is not in the graph".format(v))
        mask |= (1 << idx[v])
    return mask


def _loop_mask_from_vertices(vertices, looped_vertices):
    """
    Return bitmask of looped vertices.

    Parameters
    ----------
    vertices : list
        Ordered graph vertices.
    looped_vertices : iterable or None
        Vertices that have loops.  If None, no loops are used.
    """
    if looped_vertices is None:
        return 0
    return _bitmask_from_vertices(vertices, looped_vertices)


# ---------------------------------------------------------------------------
# Core looped zero forcing computation
# ---------------------------------------------------------------------------

def _lzf_closure(adj_mask, initial_mask, loop_mask, n):
    """
    Compute looped zero forcing closure of ``initial_mask``.

    Under the looped rule, any vertex (blue or white) can act as a forcing
    vertex if it has exactly one white neighbor, where neighbors are taken in
    the looped graph, so a looped white vertex may count itself.

    Parameters
    ----------
    adj_mask : list of int
        Underlying (simple) adjacency bitmasks.
    initial_mask : int
        Bitmask of initially blue vertices.
    loop_mask : int
        Bitmask indicating which vertices are looped.
    n : int
        Number of vertices.

    Returns
    -------
    int
        Bitmask of all blue vertices after propagation.
    """
    blue = initial_mask
    changed = True

    while changed:
        changed = False
        to_blue = 0

        for v in range(n):
            # Neighbors in the looped configuration: ordinary neighbors plus
            # itself if v is looped.
            nbrs = adj_mask[v]
            if (loop_mask >> v) & 1:
                nbrs |= (1 << v)

            white_nbrs = nbrs & ~blue

            # exactly one white neighbor
            if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                to_blue |= white_nbrs

        if to_blue:
            blue |= to_blue
            changed = True

    return blue


def _looped_zero_forcing_number_internal(adj_mask, loop_mask, n, full_mask):
    """
    Compute looped zero forcing number by increasing subset size search.

    Returns the size of the smallest looped zero forcing set found.
    """
    for size in range(0, n + 1):
        for combo in combinations(range(n), size):
            mask = 0
            for v in combo:
                mask |= (1 << v)
            if _lzf_closure(adj_mask, mask, loop_mask, n) == full_mask:
                return size
    return n  # pragma: no cover


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def looped_zero_forcing_closure(g, initial_set, looped_vertices=None):
    """
    Compute the looped zero forcing closure of ``initial_set``.

    Parameters
    ----------
    g : graph
        Graph in any format supported by this module (SageMath preferred).
    initial_set : iterable
        Vertices initially colored blue.
    looped_vertices : iterable, optional
        Vertices carrying loops in the loop configuration.

    Returns
    -------
    frozenset
        The final blue vertex set after repeated looped forcing.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    idx_to_v = {i: v for i, v in enumerate(vertices)}

    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    closure_mask = _lzf_closure(adj_mask, initial_mask, loop_mask, n)
    return frozenset(idx_to_v[i] for i in range(n) if (closure_mask >> i) & 1)


def is_looped_zero_forcing_set(g, initial_set, looped_vertices=None):
    """
    Test whether ``initial_set`` is looped zero forcing for the configuration.

    Parameters
    ----------
    g : graph
        Graph in any format supported by this module.
    initial_set : iterable
        Candidate initial blue set.
    looped_vertices : iterable, optional
        Vertices carrying loops in the loop configuration.

    Returns
    -------
    bool
        True iff ``initial_set`` forces all vertices blue.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    full_mask = (1 << n) - 1

    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    return _lzf_closure(adj_mask, initial_mask, loop_mask, n) == full_mask


def looped_zero_forcing_number(g, looped_vertices=None, return_sets=False):
    """
    Compute the looped graph zero forcing number for a fixed loop configuration.

    Parameters
    ----------
    g : graph
        A graph in any supported format (SageMath graph preferred).
    looped_vertices : iterable, optional
        Vertices carrying loops in the chosen loop configuration.
        If omitted/None, this reduces to standard zero forcing on the
        underlying graph under this module's rule implementation.
    return_sets : bool, optional
        If True, return ``(lZ, sets)`` where ``sets`` is a sorted list of
        minimum looped zero forcing sets as frozensets.

    Returns
    -------
    int
        The looped zero forcing number.
    tuple
        If ``return_sets=True``, returns ``(lZ, sets)``.
    """
    vertices, adj_mask, n = _adjacency_lists(g)

    if n == 0:
        if return_sets:
            return 0, [frozenset()]
        return 0

    full_mask = (1 << n) - 1
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    lz = _looped_zero_forcing_number_internal(adj_mask, loop_mask, n, full_mask)

    if not return_sets:
        return lz

    sets = []
    cache = {}
    for combo in combinations(range(n), lz):
        mask = 0
        for v in combo:
            mask |= (1 << v)
        if mask not in cache:
            cache[mask] = _lzf_closure(adj_mask, mask, loop_mask, n)
        if cache[mask] == full_mask:
            sets.append(frozenset(vertices[v] for v in combo))

    return lz, sorted(sets, key=lambda s: sorted(s))


def lzf(g, looped_vertices=None, return_sets=False):
    """Short alias for :func:`looped_zero_forcing_number`."""
    return looped_zero_forcing_number(
        g, looped_vertices=looped_vertices, return_sets=return_sets
    )


def load_all():
    """
    Return all public API callables exported by this module.

    Returns
    -------
    dict
        Mapping from public function name to callable.
    """
    return {
        "looped_zero_forcing_number": looped_zero_forcing_number,
        "looped_zero_forcing_closure": looped_zero_forcing_closure,
        "is_looped_zero_forcing_set": is_looped_zero_forcing_set,
        "lzf": lzf,
        "load_all": load_all,
    }
