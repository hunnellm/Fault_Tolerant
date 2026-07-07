#!/usr/bin/env python3
"""
loop_zf.py - Zero forcing utilities with explicit simple-vs-looped APIs.

Key design choice:
------------------
This module separates SIMPLE and LOOPED forcing semantics in public APIs.

- SIMPLE rule (no white forcing):
    only blue vertices may force a unique white neighbor.

- LOOPED rule (white forcing allowed):
    any vertex may force if it has exactly one white neighbor in the looped graph
    (where a specified subset of vertices has loops; this subset may be empty).

No public function infers one rule from missing parameters.
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
        raise ValueError("looped_vertices must be explicitly provided for looped functions (possibly empty set)")
    return int(_bitmask_from_vertices(vertices, looped_vertices))


# ---------------------------------------------------------------------------
# Core closures
# ---------------------------------------------------------------------------

def _zf_closure(adj_mask, initial_mask, n):
    """
    SIMPLE closure: only blue vertices can force.
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
    LOOPED closure: any vertex can force if it has a unique white neighbor
    in looped graph.
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


def _simple_zero_forcing_number_internal(adj_mask, n, full_mask):
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
# Public SIMPLE API
# ---------------------------------------------------------------------------

def simple_zero_forcing_closure(g, initial_set):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    closure_mask = _zf_closure(adj_mask, initial_mask, n)
    return frozenset(vertices[i] for i in range(n) if (closure_mask >> i) & 1)


def is_simple_zero_forcing_set(g, initial_set):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    full_mask = int((1 << n) - 1)
    return _zf_closure(adj_mask, initial_mask, n) == full_mask


def simple_zero_forcing_number(g, return_sets=False):
    vertices, adj_mask, n = _adjacency_lists(g)

    if n == 0:
        if return_sets:
            return 0, [frozenset()]
        return 0

    full_mask = int((1 << n) - 1)
    z = _simple_zero_forcing_number_internal(adj_mask, n, full_mask)
    if not return_sets:
        return z

    sets = []
    for combo in combinations(range(n), z):
        mask = 0
        for v in combo:
            mask |= (1 << int(v))
        mask = int(mask)
        if _zf_closure(adj_mask, mask, n) == full_mask:
            sets.append(frozenset(vertices[v] for v in combo))

    return z, sorted(sets, key=lambda s: sorted(s))


# ---------------------------------------------------------------------------
# Public LOOPED API
# ---------------------------------------------------------------------------

def looped_zero_forcing_closure(g, initial_set, looped_vertices):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    closure_mask = _lzf_closure(adj_mask, initial_mask, loop_mask, n)
    return frozenset(vertices[i] for i in range(n) if (closure_mask >> i) & 1)


def is_looped_zero_forcing_set(g, initial_set, looped_vertices):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    full_mask = int((1 << n) - 1)
    return _lzf_closure(adj_mask, initial_mask, loop_mask, n) == full_mask


def looped_zero_forcing_number(g, looped_vertices, return_sets=False):
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
    for combo in combinations(range(n), lz):
        mask = 0
        for v in combo:
            mask |= (1 << int(v))
        mask = int(mask)
        if _lzf_closure(adj_mask, mask, loop_mask, n) == full_mask:
            sets.append(frozenset(vertices[v] for v in combo))

    return lz, sorted(sets, key=lambda s: sorted(s))


def maximum_looped_zero_forcing_number(g, return_configurations=False, return_sets=False):
    """
    Maximum looped zero forcing number over all loop configurations.
    """
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
# Forcing record + path extraction
# ---------------------------------------------------------------------------

def _looped_forcing_record(adj_mask, initial_mask, loop_mask, n):
    n = int(n)
    full_mask = int((1 << n) - 1)
    blue = int(initial_mask)
    loop_mask = int(loop_mask)

    parent = {}
    force_edges = []

    while True:
        candidates = []  # (v,u)
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


def _simple_forcing_record(adj_mask, initial_mask, n):
    n = int(n)
    full_mask = int((1 << n) - 1)
    blue = int(initial_mask)

    parent = {}
    force_edges = []

    while True:
        candidates = []  # (v,u)
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
    Self-force u->u contributes singleton path (u,).
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

    leftovers = sorted(set(parent.keys()) - visited)
    for v in leftovers:
        if v in visited:
            continue

        if parent.get(v, None) == v:
            paths.append((vertices[v],))
            visited.add(v)
            continue

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


def simple_forcing_paths(g, initial_set, verify_minimum=False):
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, initial_set)

    if verify_minimum:
        z = simple_zero_forcing_number(g, return_sets=False)
        if bin(initial_mask).count("1") != int(z):
            raise ValueError("initial_set is not minimum size for simple zero forcing")

    full_forced, parent, _ = _simple_forcing_record(adj_mask, initial_mask, n)
    if not full_forced:
        raise ValueError("initial_set is not a simple zero forcing set")

    return _paths_from_parent(vertices, initial_mask, parent)


def looped_forcing_paths(g, initial_set, looped_vertices, verify_minimum=False):
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


# ---------------------------------------------------------------------------
# Reversal reconfiguration graphs
# Edge iff covers match path-by-path up to reversal
# ---------------------------------------------------------------------------

def _canonicalize_path_up_to_reversal(path):
    p = tuple(path)
    rp = tuple(reversed(p))
    return p if p <= rp else rp


def _path_cover_signature_strict_reversal(paths):
    canon = [_canonicalize_path_up_to_reversal(p) for p in paths]
    return tuple(sorted(canon))


def _minimum_simple_zero_forcing_sets(g):
    z, min_sets = simple_zero_forcing_number(g, return_sets=True)
    return z, min_sets


def _make_graph(vertex_labels, edges):
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


def _all_reversals_for_simple_set(g, B):
    """
    Return all reversal sets R obtainable from B over all chronological force lists.
    SIMPLE setting.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    idx = {v: i for i, v in enumerate(vertices)}
    initial_mask = _bitmask_from_vertices(vertices, B)
    full_mask = (1 << n) - 1

    if _zf_closure(adj_mask, initial_mask, n) != full_mask:
        return set()

    seqs = _all_chronological_force_lists_simple_indices(adj_mask, initial_mask, n)
    revs = set()

    for seq in seqs:
        outdeg = [0] * n
        for (u, v) in seq:
            outdeg[u] += 1
        R = frozenset(vertices[i] for i in range(n) if outdeg[i] == 0)
        revs.add(R)

    return revs


def reversal_reconfiguration_graph_simple(g, return_classes=False):
    """
    Edge B--R iff R is an achievable reversal of B from some valid chronology.
    """
    _, min_sets = _minimum_simple_zero_forcing_sets(g)
    min_sets = [frozenset(S) for S in min_sets]
    min_lookup = set(min_sets)

    rev_multi = {}  # B -> set of possible reversals among minimum sets
    for B in min_sets:
        all_rev = _all_reversals_for_simple_set(g, B)
        rev_multi[B] = {R for R in all_rev if R in min_lookup}

    edges_set = set()
    for B in min_sets:
        for R in rev_multi[B]:
            if R != B:
                e = tuple(sorted((B, R), key=lambda s: sorted(s)))
                edges_set.add(e)

    RG = _make_graph(min_sets, list(edges_set))

    if not return_classes:
        return RG

    # Connected components of this undirected reversal graph
    # (works for Sage Graph and fallback dict)
    try:
        comps = RG.connected_components()
        class_list = [sorted((frozenset(c) if isinstance(c, (set, frozenset)) else frozenset(c)) if False else list(c),
                             key=lambda s: sorted(s)) for c in comps]
    except Exception:
        # dict fallback
        seen = set()
        class_list = []
        adj = RG
        for v in sorted(adj.keys(), key=lambda s: sorted(s)):
            if v in seen:
                continue
            stack = [v]
            comp = []
            seen.add(v)
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            class_list.append(sorted(comp, key=lambda s: sorted(s)))

    class_list = sorted(class_list, key=lambda cls: [sorted(s) for s in cls])
    return RG, class_list

def _all_reversals_for_looped_set(g, B, looped_vertices):
    """
    Return all reversal sets R obtainable from B over all chronological force lists.
    LOOPED setting for fixed looped_vertices.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    initial_mask = _bitmask_from_vertices(vertices, B)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    full_mask = (1 << n) - 1

    if _lzf_closure(adj_mask, initial_mask, loop_mask, n) != full_mask:
        return set()

    seqs = _all_chronological_force_lists_looped_indices(adj_mask, initial_mask, loop_mask, n)
    revs = set()

    for seq in seqs:
        outdeg = [0] * n
        for (u, v) in seq:
            outdeg[u] += 1
        R = frozenset(vertices[i] for i in range(n) if outdeg[i] == 0)
        revs.add(R)

    return revs


def reversal_reconfiguration_graph_looped(g, looped_vertices, return_classes=False):
    """
    Edge B--R iff R is an achievable reversal of B from some valid chronology.
    LOOPED setting for fixed looped_vertices.
    """
    _, min_sets = looped_zero_forcing_number(g, looped_vertices=looped_vertices, return_sets=True)
    min_sets = [frozenset(S) for S in min_sets]
    min_lookup = set(min_sets)

    rev_multi = {}
    for B in min_sets:
        all_rev = _all_reversals_for_looped_set(g, B, looped_vertices)
        rev_multi[B] = {R for R in all_rev if R in min_lookup}

    edges_set = set()
    for B in min_sets:
        for R in rev_multi[B]:
            if R != B:
                e = tuple(sorted((B, R), key=lambda s: sorted(s)))
                edges_set.add(e)

    RG = _make_graph(min_sets, list(edges_set))

    if not return_classes:
        return RG

    try:
        comps = RG.connected_components()
        class_list = [sorted(list(c), key=lambda s: sorted(s)) for c in comps]
    except Exception:
        seen = set()
        class_list = []
        adj = RG
        for v in sorted(adj.keys(), key=lambda s: sorted(s)):
            if v in seen:
                continue
            stack = [v]
            comp = []
            seen.add(v)
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            class_list.append(sorted(comp, key=lambda s: sorted(s)))

    class_list = sorted(class_list, key=lambda cls: [sorted(s) for s in cls])
    return RG, class_list
def reversal_map_simple(g):
    """
    Return a sorted list of (B, R(B)) over minimum simple zero forcing sets.
    R(B) is computed as sinks of the forcing digraph from _simple_forcing_record.
    """
    _, min_sets = _minimum_simple_zero_forcing_sets(g)
    vertices, adj_mask, n = _adjacency_lists(g)

    out = []
    for B in min_sets:
        B = frozenset(B)
        initial_mask = _bitmask_from_vertices(vertices, B)
        full_forced, _parent, force_edges = _simple_forcing_record(adj_mask, initial_mask, n)
        if not full_forced:
            continue

        outdeg = {v: 0 for v in vertices}
        for (u_idx, _v_idx) in force_edges:
            outdeg[vertices[u_idx]] += 1

        RB = frozenset(v for v in vertices if outdeg[v] == 0)
        out.append((B, RB))

    return sorted(out, key=lambda t: (sorted(t[0]), sorted(t[1])))


def reversal_map_looped(g, looped_vertices):
    """
    Return a sorted list of (B, R(B)) over minimum looped zero forcing sets
    for the fixed loop configuration looped_vertices.
    """
    _, min_sets = looped_zero_forcing_number(g, looped_vertices=looped_vertices, return_sets=True)
    vertices, adj_mask, n = _adjacency_lists(g)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    out = []
    for B in min_sets:
        B = frozenset(B)
        initial_mask = _bitmask_from_vertices(vertices, B)
        full_forced, _parent, force_edges = _looped_forcing_record(
            adj_mask, initial_mask, loop_mask, n
        )
        if not full_forced:
            continue

        outdeg = {v: 0 for v in vertices}
        for (u_idx, _v_idx) in force_edges:
            outdeg[vertices[u_idx]] += 1

        RB = frozenset(v for v in vertices if outdeg[v] == 0)
        out.append((B, RB))

    return sorted(out, key=lambda t: (sorted(t[0]), sorted(t[1])))

# ---------------------------------------------------------------------------
# Chronological force lists (all possible)
# ---------------------------------------------------------------------------

def _single_bit_index(x):
    return int(x.bit_length() - 1)


def _possible_forces_simple(adj_mask, blue_mask, n):
    blue_mask = int(blue_mask)
    out = []
    for u in range(int(n)):
        if (blue_mask >> u) & 1:
            white_nbrs = int(adj_mask[u]) & ~blue_mask
            if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
                v = _single_bit_index(white_nbrs)
                out.append((int(u), int(v)))
    out.sort()
    return out


def _possible_forces_looped(adj_mask, blue_mask, loop_mask, n):
    blue_mask = int(blue_mask)
    loop_mask = int(loop_mask)
    out = []
    for u in range(int(n)):
        nbrs = int(adj_mask[u])
        if (loop_mask >> u) & 1:
            nbrs |= (1 << u)
        white_nbrs = nbrs & ~blue_mask
        if white_nbrs and not (white_nbrs & (white_nbrs - 1)):
            v = _single_bit_index(white_nbrs)
            out.append((int(u), int(v)))
    out.sort()
    return out


def _all_chronological_force_lists_simple_indices(adj_mask, initial_mask, n):
    n = int(n)
    initial_mask = int(initial_mask)
    full_mask = int((1 << n) - 1)
    memo = {}

    def rec(blue_mask):
        blue_mask = int(blue_mask)
        if blue_mask == full_mask:
            return [tuple()]
        if blue_mask in memo:
            return memo[blue_mask]

        candidates = _possible_forces_simple(adj_mask, blue_mask, n)
        if not candidates:
            memo[blue_mask] = []
            return memo[blue_mask]

        out = []
        for (u, v) in candidates:
            if (blue_mask >> v) & 1:
                continue
            next_mask = int(blue_mask | (1 << v))
            for tail in rec(next_mask):
                out.append(((u, v),) + tail)

        memo[blue_mask] = out
        return out

    return rec(initial_mask)


def _all_chronological_force_lists_looped_indices(adj_mask, initial_mask, loop_mask, n):
    n = int(n)
    initial_mask = int(initial_mask)
    loop_mask = int(loop_mask)
    full_mask = int((1 << n) - 1)
    memo = {}

    def rec(blue_mask):
        blue_mask = int(blue_mask)
        if blue_mask == full_mask:
            return [tuple()]
        if blue_mask in memo:
            return memo[blue_mask]

        candidates = _possible_forces_looped(adj_mask, blue_mask, loop_mask, n)
        if not candidates:
            memo[blue_mask] = []
            return memo[blue_mask]

        out = []
        for (u, v) in candidates:
            if (blue_mask >> v) & 1:
                continue
            next_mask = int(blue_mask | (1 << v))
            for tail in rec(next_mask):
                out.append(((u, v),) + tail)

        memo[blue_mask] = out
        return out

    return rec(initial_mask)


def _decorate_sequences(vertices, initial_mask, idx, seqs):
    """
    For each chronology seq=[(u1,v1),...], return triple:
      (
        seq,
        tuple(forcers_in_first_appearance_order + nonforcers),
        tuple(forced_in_order + initial_blue_vertices)
      )
    """
    results = []
    all_vertices = list(vertices)

    for seq in seqs:
        seen_forcer = set()
        forcer_order = []
        forced_order = []

        for (u, v) in seq:
            if u not in seen_forcer:
                seen_forcer.add(u)
                forcer_order.append(u)
            forced_order.append(v)

        nonforcers = [v for v in all_vertices if v not in seen_forcer]
        forcer_with_nonforcers = tuple(forcer_order + nonforcers)

        init_blues_ordered = [v for v in all_vertices if (initial_mask >> idx[v]) & 1]
        forced_plus_initial = tuple(forced_order + init_blues_ordered)

        results.append((seq, forcer_with_nonforcers, forced_plus_initial))

    return results


def all_chronological_forces_simple(g, initial_set, return_vertex_orders=False):
    vertices, adj_mask, n = _adjacency_lists(g)
    idx = {v: i for i, v in enumerate(vertices)}
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    full_mask = int((1 << n) - 1)

    if _zf_closure(adj_mask, initial_mask, n) != full_mask:
        raise ValueError("initial_set is not a simple zero forcing set")

    seqs_idx = _all_chronological_force_lists_simple_indices(adj_mask, initial_mask, n)
    seqs = [tuple((vertices[u], vertices[v]) for (u, v) in seq) for seq in seqs_idx]

    if not return_vertex_orders:
        return seqs
    return _decorate_sequences(vertices, initial_mask, idx, seqs)


def all_chronological_forces_looped(g, initial_set, looped_vertices, return_vertex_orders=False):
    vertices, adj_mask, n = _adjacency_lists(g)
    idx = {v: i for i, v in enumerate(vertices)}
    initial_mask = _bitmask_from_vertices(vertices, initial_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    full_mask = int((1 << n) - 1)

    if _lzf_closure(adj_mask, initial_mask, loop_mask, n) != full_mask:
        raise ValueError("initial_set is not a looped zero forcing set for this loop configuration")

    seqs_idx = _all_chronological_force_lists_looped_indices(adj_mask, initial_mask, loop_mask, n)
    seqs = [tuple((vertices[u], vertices[v]) for (u, v) in seq) for seq in seqs_idx]

    if not return_vertex_orders:
        return seqs
    return _decorate_sequences(vertices, initial_mask, idx, seqs)


# ---------------------------------------------------------------------------
# Reordered matrix constructors
# ---------------------------------------------------------------------------

def _base_A_matrix(g, looped_vertices=None):
    """
    Build A = adjacency(G) + 2I, then set A[v,v]=3 for v in looped_vertices.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    n = int(n)
    idx = {v: i for i, v in enumerate(vertices)}

    M = [[0] * n for _ in range(n)]
    for i in range(n):
        m = int(adj_mask[i])
        for j in range(n):
            if (m >> j) & 1:
                M[i][j] = 1

    for i in range(n):
        M[i][i] += 2

    if looped_vertices is not None:
        for v in looped_vertices:
            if v not in idx:
                raise ValueError("looped vertex {!r} is not in graph".format(v))
            M[idx[v]][idx[v]] = 3

    return vertices, idx, M


def _matrix_from_orders(vertices, idx, M, row_order, col_order):
    n = len(vertices)
    row_order = list(row_order)
    col_order = list(col_order)

    if len(row_order) != n or len(col_order) != n:
        raise ValueError("row_order and col_order must both have length n={}".format(n))
    if set(row_order) != set(vertices):
        raise ValueError("row_order is not a permutation of graph vertices")
    if set(col_order) != set(vertices):
        raise ValueError("col_order is not a permutation of graph vertices")

    row_idx = [idx[v] for v in row_order]
    col_idx = [idx[v] for v in col_order]
    R = [[M[i][j] for j in col_idx] for i in row_idx]

    try:
        return matrix(R), tuple(row_order), tuple(col_order)
    except Exception:
        return R, tuple(row_order), tuple(col_order)


def forcing_order_matrix_simple(g, force_order, include_idle_vertices=True):
    """
    Build reordered A for SIMPLE chronology.

    Rows: vertices that force in first-appearance chronological order;
          append non-forcers if include_idle_vertices=True.
    Cols: vertices forced in chronological order;
          append initial blue vertices if include_idle_vertices=True.
    """
    vertices, idx, M = _base_A_matrix(g, looped_vertices=None)
    n = len(vertices)

    force_order = list(force_order)
    for p in force_order:
        if not (isinstance(p, (tuple, list)) and len(p) == 2):
            raise ValueError("force_order must contain pairs (u,v)")
        u, v = p
        if u not in idx or v not in idx:
            raise ValueError("force pair ({!r},{!r}) uses vertex not in graph".format(u, v))

    seen_forcer = set()
    row_order = []
    col_order = []

    for (u, v) in force_order:
        if u not in seen_forcer:
            seen_forcer.add(u)
            row_order.append(u)
        col_order.append(v)

    if include_idle_vertices:
        row_order += [v for v in vertices if v not in seen_forcer]
        forced_set = set(col_order)
        initial_blues = [v for v in vertices if v not in forced_set]
        col_order += initial_blues

    if len(row_order) != n or len(col_order) != n:
        raise ValueError(
            "order lengths ({}, {}) do not match n={}; "
            "set include_idle_vertices=True for full matrix.".format(len(row_order), len(col_order), n)
        )

    return _matrix_from_orders(vertices, idx, M, row_order, col_order)


def forcing_order_matrix_looped(g, force_order, looped_vertices, include_idle_vertices=True):
    """
    Build reordered A for LOOPED chronology (loops exactly on looped_vertices).
    """
    vertices, idx, M = _base_A_matrix(g, looped_vertices=looped_vertices)
    n = len(vertices)

    force_order = list(force_order)
    for p in force_order:
        if not (isinstance(p, (tuple, list)) and len(p) == 2):
            raise ValueError("force_order must contain pairs (u,v)")
        u, v = p
        if u not in idx or v not in idx:
            raise ValueError("force pair ({!r},{!r}) uses vertex not in graph".format(u, v))

    seen_forcer = set()
    row_order = []
    col_order = []

    for (u, v) in force_order:
        if u not in seen_forcer:
            seen_forcer.add(u)
            row_order.append(u)
        col_order.append(v)

    if include_idle_vertices:
        row_order += [v for v in vertices if v not in seen_forcer]
        forced_set = set(col_order)
        initial_blues = [v for v in vertices if v not in forced_set]
        col_order += initial_blues

    if len(row_order) != n or len(col_order) != n:
        raise ValueError(
            "order lengths ({}, {}) do not match n={}; "
            "set include_idle_vertices=True for full matrix.".format(len(row_order), len(col_order), n)
        )

    return _matrix_from_orders(vertices, idx, M, row_order, col_order)


def forcing_order_matrix_from_acf_entry_simple(g, acf_entry):
    """
    Build reordered A from one entry of
      all_chronological_forces_simple(..., return_vertex_orders=True).
    """
    if not (isinstance(acf_entry, (tuple, list)) and len(acf_entry) == 3):
        raise ValueError("acf_entry must be (force_sequence, row_order, col_order)")
    _, row_order, col_order = acf_entry
    vertices, idx, M = _base_A_matrix(g, looped_vertices=None)
    return _matrix_from_orders(vertices, idx, M, row_order, col_order)


def forcing_order_matrix_from_acf_entry_looped(g, acf_entry, looped_vertices):
    """
    Build reordered A from one entry of
      all_chronological_forces_looped(..., return_vertex_orders=True).
    """
    if not (isinstance(acf_entry, (tuple, list)) and len(acf_entry) == 3):
        raise ValueError("acf_entry must be (force_sequence, row_order, col_order)")
    _, row_order, col_order = acf_entry
    vertices, idx, M = _base_A_matrix(g, looped_vertices=looped_vertices)
    return _matrix_from_orders(vertices, idx, M, row_order, col_order)


# ---------------------------------------------------------------------------
# Explicit aliases + deprecated ambiguous wrappers
# ---------------------------------------------------------------------------

def szf(g, return_sets=False):
    """Alias for simple_zero_forcing_number."""
    return simple_zero_forcing_number(g, return_sets=return_sets)


def lzf(g, looped_vertices, return_sets=False):
    """Alias for looped_zero_forcing_number."""
    return looped_zero_forcing_number(g, looped_vertices=looped_vertices, return_sets=return_sets)


def EZ(g, return_configurations=False, return_sets=False):
    """Alias for maximum_looped_zero_forcing_number."""
    return maximum_looped_zero_forcing_number(
        g,
        return_configurations=return_configurations,
        return_sets=return_sets,
    )


def zero_forcing_paths(g, initial_set, verify_minimum=False):
    """
    Deprecated ambiguous name. Use simple_forcing_paths.
    """
    return simple_forcing_paths(g, initial_set, verify_minimum=verify_minimum)


def all_chronological_forces(*args, **kwargs):
    """
    Deprecated ambiguous entry point.
    """
    raise ValueError(
        "all_chronological_forces is deprecated due to simple/looped ambiguity. "
        "Use all_chronological_forces_simple(...) or "
        "all_chronological_forces_looped(..., looped_vertices=...)."
    )


def acf(*args, **kwargs):
    """
    Deprecated ambiguous alias.
    """
    raise ValueError(
        "acf is deprecated due to simple/looped ambiguity. "
        "Use all_chronological_forces_simple(...) or "
        "all_chronological_forces_looped(..., looped_vertices=...)."
    )


def forcing_order_matrix(*args, **kwargs):
    """
    Deprecated ambiguous entry point.
    """
    raise ValueError(
        "forcing_order_matrix is deprecated due to simple/looped ambiguity. "
        "Use forcing_order_matrix_simple(...) or "
        "forcing_order_matrix_looped(..., looped_vertices=...)."
    )


def forcing_order_matrix_from_acf_entry(*args, **kwargs):
    """
    Deprecated ambiguous entry point.
    """
    raise ValueError(
        "forcing_order_matrix_from_acf_entry is deprecated due to simple/looped ambiguity. "
        "Use forcing_order_matrix_from_acf_entry_simple(...) or "
        "forcing_order_matrix_from_acf_entry_looped(..., looped_vertices=...)."
    )


def load_all():
    """
    Return public API callables.
    """
    return {
        # Simple API
        "simple_zero_forcing_closure": simple_zero_forcing_closure,
        "is_simple_zero_forcing_set": is_simple_zero_forcing_set,
        "simple_zero_forcing_number": simple_zero_forcing_number,
        "simple_forcing_paths": simple_forcing_paths,
        "all_chronological_forces_simple": all_chronological_forces_simple,
        "forcing_order_matrix_simple": forcing_order_matrix_simple,
        "forcing_order_matrix_from_acf_entry_simple": forcing_order_matrix_from_acf_entry_simple,
        "reversal_reconfiguration_graph_simple": reversal_reconfiguration_graph_simple,
        "szf": szf,

        # Looped API
        "looped_zero_forcing_closure": looped_zero_forcing_closure,
        "is_looped_zero_forcing_set": is_looped_zero_forcing_set,
        "looped_zero_forcing_number": looped_zero_forcing_number,
        "maximum_looped_zero_forcing_number": maximum_looped_zero_forcing_number,
        "looped_forcing_paths": looped_forcing_paths,
        "all_chronological_forces_looped": all_chronological_forces_looped,
        "forcing_order_matrix_looped": forcing_order_matrix_looped,
        "forcing_order_matrix_from_acf_entry_looped": forcing_order_matrix_from_acf_entry_looped,
        "reversal_reconfiguration_graph_looped": reversal_reconfiguration_graph_looped,
        "lzf": lzf,
        "EZ": EZ,

        # Deprecated wrappers / compatibility
        "zero_forcing_paths": zero_forcing_paths,
        "all_chronological_forces": all_chronological_forces,
        "acf": acf,
        "forcing_order_matrix": forcing_order_matrix,
        "forcing_order_matrix_from_acf_entry": forcing_order_matrix_from_acf_entry,

        "load_all": load_all,
    }

"""
_______________________________________________
Fort Calculations

_______________________________________________
"""
def loop_forts(g, looped_vertices, include_empty=False, include_full=True):
    """
    Return all loop forts for a fixed loop configuration.

    A loop fort S satisfies:
      for every vertex v in V(G_looped), |N(v) ∩ S| != 1,
    where G_looped has loops exactly on looped_vertices, and a loop at v means
    v is included in N(v).

    Parameters
    ----------
    g : graph-like
    looped_vertices : iterable
        Vertices that carry loops in the looped graph.
    include_empty : bool
        Whether to include the empty fort.
    include_full : bool
        Whether to include V(G) when it is a fort.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    n = int(n)

    if n == 0:
        return [frozenset()] if include_empty else []

    # Closed-neighborhood masks in the chosen loop configuration.
    nbr_masks = [0] * n
    for v in range(n):
        m = int(adj_mask[v])
        if (loop_mask >> v) & 1:
            m |= (1 << v)
        nbr_masks[v] = int(m)

    def is_loop_fort(mask):
        mask = int(mask)
        for v in range(n):
            c = int((nbr_masks[v] & mask).bit_count())
            if c == 1:
                return False
        return True

    forts = []
    full_mask = int((1 << n) - 1)

    for mask in range(1 << n):
        mask = int(mask)

        if not include_empty and mask == 0:
            continue
        if not include_full and mask == full_mask:
            continue

        if is_loop_fort(mask):
            forts.append(frozenset(vertices[i] for i in range(n) if (mask >> i) & 1))

    return sorted(forts, key=lambda s: (len(s), sorted(s)))

def is_loop_fort(g, fort_set, looped_vertices):
    """
    Decide whether fort_set is a loop fort for a fixed loop configuration.

    S is a loop fort iff for every vertex v in V(G_looped),
      |N(v) ∩ S| != 1,
    where loops are present exactly on looped_vertices.
    """
    vertices, adj_mask, n = _adjacency_lists(g)
    mask = _bitmask_from_vertices(vertices, fort_set)
    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)
    n = int(n)

    for v in range(n):
        nbrs = int(adj_mask[v])
        if (loop_mask >> v) & 1:
            nbrs |= (1 << v)
        if int((nbrs & mask).bit_count()) == 1:
            return False
    return True


def minimal_loop_forts(g, looped_vertices, include_empty=False):
    """
    Return all inclusion-minimal loop forts for a fixed loop configuration.

    By default, excludes the empty set (set include_empty=True to allow it).
    """
    forts = loop_forts(
        g,
        looped_vertices=looped_vertices,
        include_empty=include_empty,
        include_full=True,
    )

    minimal = []
    for S in forts:
        is_minimal = True
        for T in forts:
            if T != S and T.issubset(S):
                is_minimal = False
                break
        if is_minimal:
            minimal.append(S)

    return sorted(minimal, key=lambda s: (len(s), sorted(s)))

def loop_blocking_number(g, looped_vertices, return_sets=False):
    """
    Compute the minimum size of a loop blocking set for a fixed loop configuration.

    A loop blocking set S satisfies:
      for every vertex u in G (with loops added exactly on looped_vertices),
      u has at least two neighbors in S.

    Here, a loop at u contributes u itself as a neighbor when u in S.
    """
    vertices, adj_mask, n = _adjacency_lists(g)

    if n == 0:
        if return_sets:
            return 0, [frozenset()]
        return 0

    loop_mask = _loop_mask_from_vertices(vertices, looped_vertices)

    # For each vertex u, build the candidate-neighbor bitmask in the looped graph.
    nbr_masks = [0] * n
    for u in range(n):
        m = int(adj_mask[u])
        if (loop_mask >> u) & 1:
            m |= (1 << u)
        nbr_masks[u] = int(m)

    def is_loop_blocking(mask):
        mask = int(mask)
        for u in range(n):
            c = int((nbr_masks[u] & mask).bit_count())
            if c < 2:
                return False
        return True

    best_size = None
    best_sets = []

    for size in range(0, n + 1):
        found_any = False
        for combo in combinations(range(n), size):
            mask = 0
            for v in combo:
                mask |= (1 << int(v))
            mask = int(mask)

            if is_loop_blocking(mask):
                found_any = True
                if not return_sets:
                    return int(size)
                best_sets.append(frozenset(vertices[v] for v in combo))

        if found_any:
            best_size = int(size)
            break

    # If no set exists (possible for sparse/small graphs), return n and optionally [].
    if best_size is None:
        if return_sets:
            return int(n), []
        return int(n)

    return int(best_size), sorted(best_sets, key=lambda s: sorted(s))

def loop_blocking_sets(g, looped_vertices, return_sets=False):
    return loop_blocking_number(g, looped_vertices, return_sets=True)

def all_loop_configurations(g, return_vertex_sets=False):
    """
    Return all loop configurations of g as graph objects.

    Each configuration is obtained by adding loops on a subset S of V(g).
    The output is sorted deterministically by S (as sorted vertex list).

    If return_vertex_sets=True, return pairs (H, S) where:
      - H is the looped graph object
      - S is the frozenset of looped vertices used to build H
    """
    vertices, _adj_mask, n = _adjacency_lists(g)
    n = int(n)

    def _make_looped_graph(base_graph, looped_set):
        # Try Sage-style copy/add_edge(v,v)
        try:
            H = base_graph.copy()
            for v in looped_set:
                H.add_edge(v, v)
            return H
        except Exception:
            pass

        # Try NetworkX-style copy/add_edge(v,v)
        try:
            H = base_graph.copy()
            for v in looped_set:
                H.add_edge(v, v)
            return H
        except Exception:
            pass

        # Fallback: adjacency-dict graph representation with explicit loops
        verts, adj_mask_local, n_local = _adjacency_lists(base_graph)
        adj = {v: set() for v in verts}
        for i, u in enumerate(verts):
            m = int(adj_mask_local[i])
            for j, v in enumerate(verts):
                if (m >> j) & 1:
                    adj[u].add(v)
        for v in looped_set:
            adj[v].add(v)
        return {v: sorted(adj[v]) for v in sorted(adj)}

    data = []
    for mask in range(1 << n):
        mask = int(mask)
        S = frozenset(vertices[i] for i in range(n) if (mask >> i) & 1)
        H = _make_looped_graph(g, S)
        data.append((S, H))

    data = sorted(data, key=lambda t: sorted(t[0]))

    if return_vertex_sets:
        return [(H, S) for (S, H) in data]
    return [H for (_S, H) in data]

def all_loop_vertex_sets(g):
    """
    Return all loop configurations of g as vertex sets.

    Output is a deterministically ordered list of frozensets S ⊆ V(g),
    where S is interpreted as the set of looped vertices.
    """
    vertices, _adj_mask, n = _adjacency_lists(g)
    n = int(n)

    out = []
    for mask in range(1 << n):
        mask = int(mask)
        S = frozenset(vertices[i] for i in range(n) if (mask >> i) & 1)
        out.append(S)

    return sorted(out, key=lambda s: (len(s), sorted(s)))
