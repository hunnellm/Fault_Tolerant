#!/usr/bin/env python3
"""
ft_propagation.py - Fault tolerant propagation times built on top of ft_zf.py.

Definitions
-----------
Let B be a minimum 1-fault tolerant zero forcing set of a graph G.

- dfpt(G): direct fault propagation time.
  The minimum propagation time among all minimum fault tolerant zero forcing
  sets B.

- ifpt(G, B): indirect fault propagation time of a specific minimum fault
  tolerant zero forcing set B. For each v in B, compute the propagation time
  of B \ {v}; then take the maximum of those propagation times.

- ifpt(G): the minimum value of ifpt(G, B) among all minimum fault tolerant
  zero forcing sets B.

This module reuses the graph normalization and fault tolerant zero forcing
machinery from ft_zf.py.
"""

from itertools import combinations

from ft_zf import _adjacency_lists, fault_tolerant_zero_forcing_number


# ---------------------------------------------------------------------------
# Zero forcing propagation time helpers
# ---------------------------------------------------------------------------


def _bitmask_from_vertices(vertices, subset):
    """Return the bitmask for ``subset`` using the given ordered ``vertices``."""
    idx = {v: i for i, v in enumerate(vertices)}
    mask = 0
    for v in subset:
        mask |= 1 << idx[v]
    return mask



def _propagation_time_from_mask(adj_mask, initial_mask, n, full_mask):
    """
    Compute propagation time for a zero forcing set given by ``initial_mask``.

    The propagation time is the number of parallel forcing rounds needed to
    color all vertices black, where in each round every currently black vertex
    with a unique white neighbor forces that neighbor.

    Parameters
    ----------
    adj_mask : list[int]
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



def propagation_time(g, initial_set):
    """
    Compute the zero forcing propagation time of ``initial_set`` on ``g``.

    Parameters
    ----------
    g : graph
        A graph in any format supported by ``ft_zf``.
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


# ---------------------------------------------------------------------------
# Fault tolerant propagation times
# ---------------------------------------------------------------------------


def minimum_fault_tolerant_zero_forcing_sets(g, faults=1):
    """
    Return all minimum fault tolerant zero forcing sets of ``g``.

    Parameters
    ----------
    g : graph
        Graph in any format supported by ``ft_zf``.
    faults : int, optional
        Number of tolerated faults. Default is 1.

    Returns
    -------
    tuple[int, list[frozenset]]
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
        Graph in any format supported by ``ft_zf``.
    faults : int, optional
        Number of tolerated faults. Default is 1.
    return_sets : bool, optional
        If True, also return the minimum FTZF sets achieving ``dfpt``.

    Returns
    -------
    int
        The direct fault propagation time.
    tuple[int, list[frozenset]]
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
        Graph in any format supported by ``ft_zf``.
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
        Graph in any format supported by ``ft_zf``.
    faults : int, optional
        Number of tolerated faults. Default is 1.
    return_sets : bool, optional
        If True, also return the minimum FTZF sets achieving ``ifpt``.

    Returns
    -------
    int
        The indirect fault propagation time.
    tuple[int, list[frozenset]]
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


# ---------------------------------------------------------------------------
# Short aliases
# ---------------------------------------------------------------------------


def dfpt(g, faults=1, return_sets=False):
    """Short alias for :func:`direct_fault_propagation_time`."""
    return direct_fault_propagation_time(g, faults=faults, return_sets=return_sets)



def ifpt(g, faults=1, return_sets=False):
    """Short alias for :func:`indirect_fault_propagation_time`."""
    return indirect_fault_propagation_time(g, faults=faults, return_sets=return_sets)
