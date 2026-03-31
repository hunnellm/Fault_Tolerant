# Fault Tolerant Zero Forcing

Python implementation of fault tolerant zero forcing numbers for graphs.

## Overview

**Zero forcing** is a graph colouring process: start with a set *S* of initially coloured (black) vertices; repeatedly apply the rule *"a black vertex v forces its unique white neighbour u to become black"* until no more forces can be made.  The **zero forcing number** Z(G) is the minimum size of a set *S* that eventually colours every vertex.

A **k-fault tolerant zero forcing set** relaxes the requirement so that up to *k* vertices in *S* may fail: *S* is a k-FTZF set if every subset of *S* with |S|−k elements is a zero forcing set.  The **fault tolerant zero forcing number** ftZ(G, k) is the minimum size of any k-FTZF set.

When k = 0, ftZ(G, 0) = Z(G).

---

## Public API (`ft_zf.py`)

### `fault_tolerant_zero_forcing_number(g, faults=1, return_sets=False)`

Compute the fault tolerant zero forcing number ftZ(G, faults).

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `g` | graph | — | Input graph (see *Graph formats* below). |
| `faults` | `int` | `1` | Number of faults to tolerate (k). Must be ≥ 0. |
| `return_sets` | `bool` | `False` | If `True`, return all minimum k-fault tolerant zero forcing sets alongside the number. |

**Returns**

- `int` — the fault tolerant zero forcing number (when `return_sets=False`).
- `(int, list[frozenset])` — a tuple `(ftZ, sets)` where *sets* is a sorted list of `frozenset` objects, one per minimum k-FTZF set (when `return_sets=True`).
- `-1` — when no k-FTZF set exists (e.g. k ≥ n for a non-empty graph).

**Raises**

- `ValueError` — if `faults < 0`.

---

### `zero_forcing_number(g)`

Compute the standard zero forcing number Z(G).  Equivalent to
`fault_tolerant_zero_forcing_number(g, faults=0)`.

---

### `ftz(g, faults=1, return_sets=False)`

Short alias for `fault_tolerant_zero_forcing_number`.

---

### `Z(g)`

Short alias for `zero_forcing_number`.

---

### `load_all()`

Return a `dict` mapping every public name in this module to its callable.
Calling `load_all()` is safe at any time and is fully idempotent — repeated
calls always return an equivalent mapping.

**Returns**

| Key | Value |
|---|---|
| `"fault_tolerant_zero_forcing_number"` | compute ftZ(G, k) |
| `"zero_forcing_number"` | compute Z(G) |
| `"ftz"` | alias for `fault_tolerant_zero_forcing_number` |
| `"Z"` | alias for `zero_forcing_number` |
| `"load_all"` | this function |

```python
from ft_zf import load_all
import networkx as nx

api = load_all()
api["zero_forcing_number"](nx.path_graph(5))           # 1
api["fault_tolerant_zero_forcing_number"](nx.path_graph(5))  # 2
api["ftz"](nx.path_graph(5), faults=0)                 # 1
```

---

### Graph formats

All functions accept three graph representations:

| Format | Example |
|---|---|
| **NetworkX** `Graph` | `networkx.path_graph(5)` |
| **SageMath** graph | `graphs.PathGraph(5)` |
| **dict** `{v: [neighbours]}` | `{0: [1], 1: [0, 2], 2: [1]}` |

---

## Example usage

```python
import networkx as nx
from ft_zf import fault_tolerant_zero_forcing_number, zero_forcing_number

# Path graph P_5 = 0-1-2-3-4
g = nx.path_graph(5)

# Standard zero forcing number
print(zero_forcing_number(g))          # 1

# 1-fault tolerant zero forcing number (default)
print(fault_tolerant_zero_forcing_number(g))           # 2

# Same as faults=0 → Z(G)
print(fault_tolerant_zero_forcing_number(g, faults=0)) # 1

# Return all minimum 1-FT ZF sets
num, sets = fault_tolerant_zero_forcing_number(g, faults=1, return_sets=True)
print(num)                             # 2
print(sets)                            # [frozenset({0, 4})]

# 2-fault tolerant zero forcing number
print(fault_tolerant_zero_forcing_number(g, faults=2)) # 4

# Plain dict input
g_dict = {0: [1], 1: [0, 2], 2: [1]}  # P_3
print(fault_tolerant_zero_forcing_number(g_dict))      # 2
```

### Known values

| Graph | Z(G) | ftZ(G, 1) |
|---|---|---|
| Path P_n (n ≥ 2) | 1 | 2 |
| Cycle C_4 | 2 | 4 |
| Cycle C_5 | 2 | 4 |
| Complete K_n | n−1 | n |

---

## Algorithm & complexity

1. **Zero forcing closure** — computed with a bitmask representation (each vertex set is an integer); the colour-change rule is applied by checking whether a single bit is set in the white-neighbour mask.  Time: O(n²) per closure.

2. **Zero forcing number** — found by iterating over subsets of increasing size.

3. **Fault tolerant search** — uses the lower bound ftZ(G, k) ≥ Z(G) + k to prune the search space, starting at size Z(G) + k.  For each candidate set S, all C(|S|, k) sub-sets of size |S|−k are checked as zero forcing sets.

4. **Memoisation** — a shared cache maps each integer bitmask to its zero forcing closure.  This is shared across all subset checks (inner and outer), so repeated work is avoided.

5. **Early exit** — when `return_sets=False`, the search terminates as soon as the first minimum set is found.

**Worst-case complexity**: O(C(n, ftZ) · C(ftZ, k) · n²) with the memoisation cache reducing the effective number of zero forcing closure calls.

**Practical limits**: graphs with up to ~20 vertices complete in seconds for faults ≤ 2.  Larger graphs are feasible when ftZ is small relative to n.  The bitmask representation handles up to ~60 vertices correctly (Python integers have arbitrary precision).

---

## Running tests

```bash
pip install networkx pytest
pytest test_ft_zf.py -v
```

---

## File overview

| File | Description |
|---|---|
| `ft_zf.py` | Fast pure-Python implementation (no SageMath required). |
| `test_ft_zf.py` | pytest test suite covering correctness and edge cases. |
| `fault_tolerant.py` | Original SageMath-based zero forcing functions. |
| `helper.py` | SageMath helper functions (loads external Sage library). |
| `misc_fncs.py` | Miscellaneous utility functions. |
