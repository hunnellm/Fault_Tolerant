#!/usr/bin/env python3
"""
loop_zf.py - Looped and regular zero forcing utilities.

This module provides:
- Looped zero forcing number for a fixed loop configuration.
- Maximum looped zero forcing number over all loop configurations.
- Forcing-path extraction (regular and looped).
- Reversal reconfiguration graphs (regular and looped), where two minimum
  forcing sets are adjacent iff their path covers are pairwise reversals:
    each path in one cover must match a reversed path in the other cover
    (equivalently, they match after canonicalizing each path up to reversal).

All bitmask arithmetic is coerced to Python int for Sage safety.
"""

from itertools import combinations


# ---------------------------------------------------------------------------
# Graph normalization
# ---------------------------------------------------------------------------

def _adjacency_lists(g):
    """
    Return (vertices, adj_mask, n) in bitmask-ready format.
    """
    if isinstance(g, dict):
        vertices = sorted(g.keys())
        raw_adj = {v: list(g[v]) for v in vertices}
    elif hasattr(g, "adjacency"):
        # NetworkX-like
        adj_raw = dict(g.adjacency())
        vertices = sorted(adj_raw.keys())
        raw_adj = {v: list(adj_raw[v].keys()) for v in vertices}
    else:
        # SageMath-compatible
        vertices = sorted(g.vertices())
        raw_adj = {v: list(g.neighbors(v)) for v in vertices}

    n = int(len(vertices))
    idx = {v: int(i) for i, v in enumerate(vertices)}

    adj_mask = [0] * n
    for i, v in enumerate(vertices):
        m = 0
        for u in raw_adj[v]:
            if u in idx:
                m |= (1 << idx[u])
        adj_mask[int(i)] = int(m)

    return vertices, adj_mask, n


def _bitmask_from_vertices(vertices, subset):
    idx = {v: int(i) for i, v in enumerate(vertices)}
    mask = 0
    for v in subset:
        if v not in idx:
            raise ValueError("vertex {!r} is not in the graph".format(v))
        mask |= (1 << idx[v])
    return int(mask)


def _loop_mask_from_vertices(vertices, looped_vertices):
    if looped_vertices is None:
        return 0
    return int(_bitmask_from_vertices(vertices, looped_vertices))


# ---------------------------------------------------------------------------
# Core closures
# ---------------------------------------------------------------------------

def _zf_closure(adj_mask, initial_mask, n):
    """
    Regular zero forcing closure.
    Only blue vertices may force their unique white neighbor.
    """
    black = int(initial_mask)
    n = int(n)

    changed = True
    while changed:
        changed = False
        for v in range(n):
            if (black >> v) & 1:
                white_nbrs = int(adj_mask[v]) & ~black
                if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                    black |= int(white_nbrs)
                    changed = True
    return int(black)


def _lzf_closure(adj_mask, initial_mask, loop_mask, n):
    """
    Looped zero forcing closure.
    Any vertex may force if it has exactly one white neighbor in looped graph.
    """
    blue = int(initial_mask)
    loop_mask = int(loop_mask)
    n = int(n)

    changed = True
    while changed:
        changed = False
        to_blue = 0

        for v in range(n):
            nbrs = int(adj_mask[v])
            if (loop_mask >> v) & 1:
                nbrs |= (1 << v)

            white_nbrs = nbrs & ~blue
            if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                to_blue |= int(white_nbrs)

        if to_blue:
            blue |= int(to_blue)
            changed = True

    return int(blue)


def _zero_forcing_number_internal(adj_mask, n, full_mask):
    n = int(n)
    full_mask = int(full_mask)
    for size in range(0, n + 1):
        for combo in combinations(range(n), size):
            mask = 0
            for v in combo:
                mask |= (1 << int(v))
            mask = int(mask)
            if _zf_closure(adj_mask, mask, n) == full_mask:
                return int(size)
    return int(n)


def _looped_zero_forcing_number_internal(adj_mask, loop_mask, n, full_mask):
    n = int(n)
    loop_mask = int(loop_mask)
    full_mask = int(full_mask)

    for size in range(0, n + 1):
        for combo in combinations(range(n), size):
            mask = 0
            for v in combo:
                mask |= (1 << int(v))
            mask = int(mask)
            if _lzf_closure(adj_mask, mask, loop_mask, n) == full_mask:
                return int(size)
    return int(n)


# ---------------------------------------------------------------------------
# Public: fixed loop configuration
# ---------------------------------------------------------------------------

def looped_zero_forcing_closure(g, initial_set, looped_vertices=None):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    closure_mask = _lzf_closure(adj_mask, initial_mask, loop_mask, n)
    return frozenset(vertices[i] for i in range(n) if (closure_mask >> i) & 1)


def is_looped_zero_forcing_set(g, initial_set, looped_vertices=None):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    full_mask = int((1 << n) - 1)
    return _lzf_closure(adj_mask, initial_mask, loop_mask, n) == full_mask


def looped_zero_forcing_number(g, looped_vertices=None, return_sets=False):
    vertices, adj_mask, n = _adjacency_lists(g)

    if n == 0:
        if return_sets:
            return 0, [frozenset()]
        return 0

    full_mask = int((1 << n) - 1)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    lz = _looped_zero_forcing_number_internal(adj_mask, loop_mask, n, full_mask)
    if not return_sets:
        return lz

    sets = []
    cache = {}
    for combo in combinations(range(n), lz):
        mask = 0
        for v in combo:
            mask |= (1 << int(v))
        mask = int(mask)
        if mask not in cache:
            cache[mask] = _lzf_closure(adj_mask, mask, loop_mask, n)
        if cache[mask] == full_mask:
            sets.append(frozenset(vertices[v] for v in combo))

    return lz, sorted(sets, key=lambda s: sorted(s))


# ---------------------------------------------------------------------------
# Public: max over all loop configurations
# ---------------------------------------------------------------------------

def maximum_looped_zero_forcing_number(g, return_configurations=False, return_sets=False):
    vertices, adj_mask, n = _adjacency_lists(g)

    if n == 0:
        if return_sets:
            return 0, [(frozenset(), [frozenset()])]
        if return_configurations:
            return 0, [frozenset()]
        return 0

    full_mask = int((1 << n) - 1)
    max_lz = -1
    maximizing_configs = []

    for loop_mask in range(1 << int(n)):
        loop_mask = int(loop_mask)
        lz = _looped_zero_forcing_number_internal(adj_mask, loop_mask, n, full_mask)

        if lz > max_lz:
            max_lz = lz
            maximizing_configs = [loop_mask]
        elif lz == max_lz:
            maximizing_configs.append(loop_mask)

    def mask_to_set(m):
        return frozenset(vertices[i] for i in range(n) if (m >> i) & 1)

    if not return_configurations and not return_sets:
        return max_lz

    if return_configurations and not return_sets:
        cfgs = [mask_to_set(m) for m in maximizing_configs]
        return max_lz, sorted(cfgs, key=lambda s: sorted(s))

    data = []
    for m in maximizing_configs:
        cfg = mask_to_set(m)
        _, min_sets = looped_zero_forcing_number(g, looped_vertices=cfg, return_sets=True)
        data.append((cfg, min_sets))

    data = sorted(data, key=lambda pair: sorted(pair[0]))
    return max_lz, data


# ---------------------------------------------------------------------------
# Forcing-path extraction (looped and regular)
# ---------------------------------------------------------------------------

def _looped_forcing_record(adj_mask, initial_mask, loop_mask, n):
    """
    Record looped forcing map parent[v] = u where u forced v.

    Tie-breaking in same round for determinism:
      smallest forced vertex index v first, then smallest u.
    """
    n = int(n)
    full_mask = int((1 << n) - 1)
    blue = int(initial_mask)
    loop_mask = int(loop_mask)

    parent = {}
    force_edges = []

    while True:
        candidates = []  # (v, u)
        for u in range(n):
            nbrs = int(adj_mask[u])
            if (loop_mask >> u) & 1:
                nbrs |= (1 << u)

            white_nbrs = nbrs & ~blue
            if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                v = int(white_nbrs.bit_length() - 1)
                candidates.append((v, u))

        if not candidates:
            break

        candidates.sort()
        to_blue = 0
        used_v = set()
        round_edges = []

        for v, u in candidates:
            if v in used_v:
                continue
            used_v.add(v)
            to_blue |= (1 << v)
            round_edges.append((u, v))

        if to_blue == 0:
            break

        for (u, v) in round_edges:
            if v not in parent:
                parent[v] = u
                force_edges.append((u, v))

        blue |= int(to_blue)

    return bool(blue == full_mask), parent, force_edges


def _zf_forcing_record(adj_mask, initial_mask, n):
    """
    Record regular forcing map parent[v] = u where u forced v.
    """
    n = int(n)
    full_mask = int((1 << n) - 1)
    blue = int(initial_mask)

    parent = {}
    force_edges = []

    while True:
        candidates = []  # (v, u)
        for u in range(n):
            if (blue >> u) & 1:
                white_nbrs = int(adj_mask[u]) & ~blue
                if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                    v = int(white_nbrs.bit_length() - 1)
                    candidates.append((v, u))

        if not candidates:
            break

        candidates.sort()
        to_blue = 0
        used_v = set()
        round_edges = []

        for v, u in candidates:
            if v in used_v:
                continue
            used_v.add(v)
            to_blue |= (1 << v)
            round_edges.append((u, v))

        if to_blue == 0:
            break

        for (u, v) in round_edges:
            if v not in parent:
                parent[v] = u
                force_edges.append((u, v))

        blue |= int(to_blue)

    return bool(blue == full_mask), parent, force_edges


def _paths_from_parent(vertices, initial_mask, parent):
    """
    Build directed forcing paths from parent map.
    Self-force u->u becomes singleton path (u,).
    """
    n = len(vertices)

    children = {}
    for v, u in parent.items():
        if u == v:
            children.setdefault(v, [])
            continue
        children.setdefault(u, []).append(v)

    for u in children:
        children[u].sort()

    starts = [i for i in range(n) if (initial_mask >> i) & 1]
    starts.sort()

    visited = set()
    paths = []

    # Build paths from initial blue vertices first
    for s in starts:
        if s in visited:
            continue

        if parent.get(s, None) == s:
            paths.append((vertices[s],))
            visited.add(s)
            continue

        path = [s]
        visited.add(s)
        cur = s
        while True:
            nxts = children.get(cur, [])
            if len(nxts) != 1:
                break
            nxt = nxts[0]
            if nxt in visited:
                break
            path.append(nxt)
            visited.add(nxt)
            cur = nxt

        paths.append(tuple(vertices[i] for i in path))

    # Any uncovered forced vertices
    leftovers = sorted(set(parent.keys()) - visited)
    for v in leftovers:
        if v in visited:
            continue

        if parent.get(v, None) == v:
            paths.append((vertices[v],))
            visited.add(v)
            continue

        # Backtrack to root
        root = v
        seen = set()
        while root in parent and parent[root] != root and root not in seen:
            seen.add(root)
            root = parent[root]

        if root in visited:
            continue

        path = [root]
        visited.add(root)
        cur = root
        while True:
            nxts = children.get(cur, [])
            if len(nxts) != 1:
                break
            nxt = nxts[0]
            if nxt in visited:
                break
            path.append(nxt)
            visited.add(nxt)
            cur = nxt

        paths.append(tuple(vertices[i] for i in path))

    return sorted(paths, key=lambda p: (len(p), p))


def looped_forcing_paths(g, initial_set, looped_vertices=None, verify_minimum=False):
    """
    Return directed forcing paths for looped forcing.
    If a white vertex forces itself, path is singleton (v,).
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    if verify_minimum:
        lz = looped_zero_forcing_number(g, looped_vertices=looped_vertices, return_sets=False)
        if bin(initial_mask).count("1") != int(lz):
            raise ValueError("initial_set is not minimum size for this loop configuration")

    full_forced, parent, _ = _looped_forcing_record(adj_mask, initial_mask, loop_mask, n)
    if not full_forced:
        raise ValueError("initial_set is not a looped zero forcing set")

    return _paths_from_parent(vertices, initial_mask, parent)


def zero_forcing_paths(g, initial_set, verify_minimum=False):
    """
    Return directed forcing paths for regular zero forcing.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)

    if verify_minimum:
        z = _zero_forcing_number_internal(adj_mask, n, int((1 << n) - 1))
        if bin(initial_mask).count("1") != int(z):
            raise ValueError("initial_set is not minimum size for regular zero forcing")

    full_forced, parent, _ = _zf_forcing_record(adj_mask, initial_mask, n)
    if not full_forced:
        raise ValueError("initial_set is not a zero forcing set")

    return _paths_from_parent(vertices, initial_mask, parent)


# ---------------------------------------------------------------------------
# Reversal logic: strict pairwise path reversal
# ---------------------------------------------------------------------------

def _canonicalize_path_up_to_reversal(path):
    """
    Canonical representative of a path tuple up to reversal.
    """
    p = tuple(path)
    rp = tuple(reversed(p))
    return p if p <= rp else rp


def _path_cover_signature_strict_reversal(paths):
    """
    Canonical signature for path cover under strict pairwise reversal condition.

    Two covers have the same signature iff each path in one corresponds to the
    reverse of a path in the other (allowing identical singleton paths).
    This is equivalent to multiset equality after canonicalizing each path up to
    reversal and sorting.
    """
    canon = [_canonicalize_path_up_to_reversal(p) for p in paths]
    return tuple(sorted(canon))


# ---------------------------------------------------------------------------
# Reversal reconfiguration graphs
# ---------------------------------------------------------------------------

def _minimum_zero_forcing_sets(g):
    vertices, adj_mask, n = _adjacency_lists(g)
    if n == 0:
        return 0, [frozenset()]

    full_mask = int((1 << n) - 1)
    z = _zero_forcing_number_internal(adj_mask, n, full_mask)

    sets = []
    for combo in combinations(range(n), z):
        mask = 0
        for v in combo:
            mask |= (1 << int(v))
        mask = int(mask)

        if _zf_closure(adj_mask, mask, n) == full_mask:
            sets.append(frozenset(vertices[v] for v in combo))

    return z, sorted(sets, key=lambda s: sorted(s))


def _make_graph(vertex_labels, edges):
    """
    Return Sage Graph when available; fallback adjacency dict otherwise.
    """
    try:
        H = Graph()
        H.add_vertices(vertex_labels)
        H.add_edges(edges)
        return H
    except Exception:
        adj = {v: set() for v in vertex_labels}
        for u, v in edges:
            adj[u].add(v)
            adj[v].add(u)
        return {v: sorted(adj[v]) for v in sorted(adj)}


def reversal_reconfiguration_graph_simple(g, return_classes=False):
    """
    Reversal reconfiguration graph for regular zero forcing.

    Vertices:
        minimum zero forcing sets of G.
    Edges:
        between two vertices iff every path in one path cover is the reverse
        of a path in the other (pairwise, as multisets).
    """
    _, min_sets = _minimum_zero_forcing_sets(g)

    classes = {}
    for S in min_sets:
        sig = _path_cover_signature_strict_reversal(
            zero_forcing_paths(g, S, verify_minimum=False)
        )
        classes.setdefault(sig, []).append(S)

    verts = list(min_sets)
    edges = []
    for group in classes.values():
        if len(group) >= 2:
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    edges.append((group[i], group[j]))

    RG = _make_graph(verts, edges)

    if return_classes:
        class_list = [sorted(v, key=lambda s: sorted(s)) for v in classes.values()]
        class_list = sorted(class_list, key=lambda cls: [sorted(s) for s in cls])
        return RG, class_list
    return RG


def reversal_reconfiguration_graph_looped(g, looped_vertices=None, return_classes=False):
    """
    Reversal reconfiguration graph for looped forcing (fixed loop configuration).

    Vertices:
        minimum looped zero forcing sets for the given loop configuration.
    Edges:
        between two vertices iff every path in one looped path cover is the
        reverse of a path in the other (pairwise, as multisets).
    """
    _, min_sets = looped_zero_forcing_number(
        g, looped_vertices=looped_vertices, return_sets=True
    )

    classes = {}
    for S in min_sets:
        sig = _path_cover_signature_strict_reversal(
            looped_forcing_paths(
                g, S, looped_vertices=looped_vertices, verify_minimum=False
            )
        )
        classes.setdefault(sig, []).append(S)

    verts = list(min_sets)
    edges = []
    for group in classes.values():
        if len(group) >= 2:
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    edges.append((group[i], group[j]))

    RG = _make_graph(verts, edges)

    if return_classes:
        class_list = [sorted(v, key=lambda s: sorted(s)) for v in classes.values()]
        class_list = sorted(class_list, key=lambda cls: [sorted(s) for s in cls])
        return RG, class_list
    return RG


# ---------------------------------------------------------------------------
# Aliases + exports
# ---------------------------------------------------------------------------

def lzf(g, looped_vertices=None, return_sets=False):
    """Alias for looped_zero_forcing_number."""
    return looped_zero_forcing_number(
        g, looped_vertices=looped_vertices, return_sets=return_sets
    )


def EZ(g, return_configurations=False, return_sets=False):
    """Alias for maximum_looped_zero_forcing_number."""
    return maximum_looped_zero_forcing_number(
        g,
        return_configurations=return_configurations,
        return_sets=return_sets,
    )


def load_all():
    """
    Return public API callables.
    """
    return {
        "looped_zero_forcing_number": looped_zero_forcing_number,
        "looped_zero_forcing_closure": looped_zero_forcing_closure,
        "is_looped_zero_forcing_set": is_looped_zero_forcing_set,
        "maximum_looped_zero_forcing_number": maximum_looped_zero_forcing_number,
        "looped_forcing_paths": looped_forcing_paths,
        "zero_forcing_paths": zero_forcing_paths,
        "reversal_reconfiguration_graph_simple": reversal_reconfiguration_graph_simple,
        "reversal_reconfiguration_graph_looped": reversal_reconfiguration_graph_looped,
        "lzf": lzf,
        "EZ": EZ,
        "load_all": load_all,
    }
