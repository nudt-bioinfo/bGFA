"""
Generate simple GFA graphs for testing.

This script creates DNA sequences and graph edges with configurable
complexity, optionally compressing linear chains into single segments.

Author: grampart
Version: 0.1
Copyright: Copyright(c) 2025 by grampart
License: MIT License
"""
# Generate a simple GFA file: configurable node count, sequence length, and graph complexity

import argparse
import os
import random
import sys
from collections import defaultdict, deque
from typing import Dict, List, Tuple

DNA_ALPHABET = "ACGT"


def random_dna(length: int, rng: random.Random) -> str:
    """Generate a random DNA string of 90~110% length using provided RNG."""
    # Choose an actual length uniformly in [0.9*length, 1.1*length]
    min_len = max(1, int(length * 0.9))
    max_len = max(min_len, int(length * 1.1))
    real_len = rng.randint(min_len, max_len)
    return "".join(rng.choice(DNA_ALPHABET) for _ in range(real_len))


def make_edges(
    num_nodes: int,
    complexity: int,
    loop_prob: float,
    rng: random.Random,
) -> List[Tuple[int, int, str, str]]:
    """Construct edges for a test graph.

    Returns a list of edges in the form:
    (from_id, to_id, from_orient, to_orient), with orientations in '+' or '-'.
    """

    edges: List[Tuple[int, int, str, str]] = []

    # Base linear backbone
    for node_id in range(1, num_nodes):
        edges.append((node_id, node_id + 1, "+", "+"))

    # Add branches and loops per complexity
    extra_edges = 0
    if complexity <= 0:
        return edges
    elif complexity == 1:
        extra_edges = max(1, num_nodes // 4)
    elif complexity == 2:
        extra_edges = max(2, num_nodes // 3)
    else:
        extra_edges = max(3, num_nodes // 2)

    # Helper: check if adding edge (a->b) will create a cycle (loop)
    def would_create_cycle(src: int, dst: int, current_edges: List[Tuple[int, int, str, str]]) -> bool:
        """Return True if adding edge src->dst forms a cycle.

        A cycle exists when src == dst or there is already a path dst->...->src.
        """
        if src == dst:
            return True
        # BFS/DFS to see if there is a path from b to a in current graph
        adj = defaultdict(list)
        for u, v, _oa, _ob in current_edges:
            adj[u].append(v)
        seen = set()
        dq = deque([dst])
        while dq:
            v = dq.popleft()
            if v == src:
                return True
            if v in seen:
                continue
            seen.add(v)
            for w in adj.get(v, []):
                if w not in seen:
                    dq.append(w)
        return False

    for _ in range(extra_edges):
        a = rng.randint(1, num_nodes)
        b = rng.randint(1, num_nodes)
        # Decide if we add a cycle-forming edge
        if would_create_cycle(a, b, edges):
            if rng.random() < loop_prob:
                edges.append(
                    (a, b, rng.choice(["+", "-"]), rng.choice(["+", "-"]))
                )
            # else: skip cycle edge this time
        else:
            # Non-cycle edge: add directly (no extra branch_prob gating)
            edges.append(
                (a, b, rng.choice(["+", "-"]), rng.choice(["+", "-"])))

    return edges


def _is_strict_linear(num_nodes: int, edges: List[Tuple[int, int, str, str]]) -> bool:
    """Check if the whole graph forms a strict linear chain."""
    if num_nodes == 0:
        return False
    # Must have exactly n-1 edges
    if len(edges) != max(0, num_nodes - 1):
        return False
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    nodes = set(range(1, num_nodes + 1))
    for a, b, _oa, _ob in edges:
        outdeg[a] += 1
        indeg[b] += 1
    # Degree profile checks
    starts = [v for v in nodes if indeg[v] == 0 and outdeg[v] == 1]
    ends = [v for v in nodes if indeg[v] == 1 and outdeg[v] == 0]
    mids_ok = all((indeg[v] == 1 and outdeg[v] == 1) or (
        indeg[v] + outdeg[v] == 0) for v in nodes if v not in set(starts + ends))
    if len(starts) != 1 or len(ends) != 1 or not mids_ok:
        return False
    # Connectivity check
    s = starts[0]
    seen = set([s])
    cur = s
    for _ in range(num_nodes - 1):
        nxt = None
        for a, b, _oa, _ob in edges:
            if a == cur:
                nxt = b
                break
        if nxt is None or nxt in seen:
            return False
        seen.add(nxt)
        cur = nxt
    return len(seen) == num_nodes


def write_gfa(output_path: str, sequences: List[str], edges: List[Tuple[int, int, str, str]]) -> None:
    """Write GFA v1: S/L records; add one P record if strictly linear.

    If `output_path` is '-', emit to stdout instead of a file.
    """
    if output_path == "-":
        f = sys.stdout
        # Avoid mixing binary and text modes; stdout is fine
        f.write("H\tVN:Z:1\n")
        for i, seq in enumerate(sequences, start=1):
            f.write(f"S\t{i}\t{seq}\n")
        for (a, b, oa, ob) in edges:
            f.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M\n")
        if len(sequences) > 0 and _is_strict_linear(len(sequences), edges):
            path_nodes = ",".join(str(i) for i in range(1, len(sequences) + 1))
            f.write(f"P\tlinear\t{path_nodes}\t*\n")
        return

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("H\tVN:Z:1\n")
        for i, seq in enumerate(sequences, start=1):
            f.write(f"S\t{i}\t{seq}\n")
        for (a, b, oa, ob) in edges:
            f.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M\n")
        if len(sequences) > 0 and _is_strict_linear(len(sequences), edges):
            path_nodes = ",".join(str(i) for i in range(1, len(sequences) + 1))
            f.write(f"P\tlinear\t{path_nodes}\t*\n")


def compress_linear(
    sequences: List[str],
    edges: List[Tuple[int, int, str, str]],
) -> Tuple[List[str], List[Tuple[int, int, str, str]], Dict[int, int]]:
    """
    Merge maximal strictly linear chains where internal nodes have in-degree==1 and out-degree==1.
    Returns:
    - New sequence list (reindexed 1..M)
    - New edge list (orientations normalized to '+')
    - Mapping from old id to new id
    Note: For simplicity and robustness, orientations are normalized to '+', and overlap stays as '0M'.
    """
    num_nodes = len(sequences)
    id_to_seq: Dict[int, str] = {i + 1: s for i, s in enumerate(sequences)}
    outgoing: Dict[int, List[Tuple[int, int, str, str]]] = defaultdict(list)
    incoming: Dict[int, List[Tuple[int, int, str, str]]] = defaultdict(list)
    for e in edges:
        a, b, oa, ob = e
        outgoing[a].append(e)
        incoming[b].append(e)

    indeg = {i: len(incoming.get(i, [])) for i in range(1, num_nodes + 1)}
    outdeg = {i: len(outgoing.get(i, [])) for i in range(1, num_nodes + 1)}

    visited = set()
    chains: List[List[int]] = []

    # Find chain starts: out==1 and in!=1, or sources with in==0
    for v in range(1, num_nodes + 1):
        if v in visited and indeg.get(v, 0) == 0 and outdeg.get(v, 0) == 0:
            continue
        if outdeg.get(v, 0) == 1 and indeg.get(v, 0) != 1:
            cur = v
            chain = [cur]
            visited.add(cur)
            while outdeg.get(cur, 0) == 1:
                e = outgoing[cur][0]
                nxt = e[1]
                if indeg.get(nxt, 0) != 1 or nxt in visited:
                    break
                chain.append(nxt)
                visited.add(nxt)
                cur = nxt
            chains.append(chain)

    # Handle unvisited segments with out==1 and in==1 (find leftmost and expand)
    for v in range(1, num_nodes + 1):
        if v in visited:
            continue
        if outdeg.get(v, 0) == 1 and indeg.get(v, 0) == 1:
            # Backtrack to the leftmost
            left = v
            while indeg.get(left, 0) == 1:
                prev_e = incoming[left][0]
                prev_v = prev_e[0]
                if outdeg.get(prev_v, 0) != 1 or prev_v in visited:
                    break
                left = prev_v
            # Expand to the right from left
            cur = left
            chain = [cur]
            visited.add(cur)
            while outdeg.get(cur, 0) == 1:
                e = outgoing[cur][0]
                nxt = e[1]
                if indeg.get(nxt, 0) != 1 or nxt in visited:
                    break
                chain.append(nxt)
                visited.add(nxt)
                cur = nxt
            chains.append(chain)

    # Remaining isolated nodes (out==0 and in==0) are also chains
    for v in range(1, num_nodes + 1):
        if v not in visited:
            chains.append([v])
            visited.add(v)

    # Build new graph
    new_id = 0
    old_to_new: Dict[int, int] = {}
    new_seq_by_id: Dict[int, str] = {}
    for chain in chains:
        new_id += 1
        merged_seq = "".join(id_to_seq[x] for x in chain)
        new_seq_by_id[new_id] = merged_seq
        for x in chain:
            old_to_new[x] = new_id

    # Rebuild edges: keep inter-chain edges; remove internal; normalize orient to '+'
    new_edges_set = set()
    for a, b, _oa, _ob in edges:
        na, nb = old_to_new[a], old_to_new[b]
        if na == nb:
            continue
        new_edges_set.add((na, nb, "+", "+"))
    new_edges = list(new_edges_set)

    # Sort by new id and emit lists
    sequences_out = [new_seq_by_id[i] for i in range(1, new_id + 1)]
    return sequences_out, new_edges, old_to_new


def parse_args(argv: List[str]) -> argparse.Namespace:
    """Parse CLI arguments for graph generation."""
    p = argparse.ArgumentParser(description="Generate a simple GFA graph file")
    p.add_argument("-n", "--nodes", type=int,
                   default=10, help="Number of nodes")
    p.add_argument("-l", "--length", type=int, default=100,
                   help="Sequence length of each node")
    p.add_argument("-c", "--complexity", type=int, default=1,
                   choices=[0, 1, 2, 3], help="Graph complexity level 0-3")
    p.add_argument("-o", "--output", type=str,
                   default="temp/generated.gfa", help="Output GFA file path")
    p.add_argument("--emit", type=str, choices=["file", "stdout"], default="file",
                   help="Emit target: save to file or print to stdout ('-' path also prints)")
    p.add_argument("--loop-prob", type=float,
                   default=0.2, help="Probability of loops (a->a)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed (reproducible)")
    return p.parse_args(argv)


def main(argv: List[str]) -> int:
    """Entry point: build sequences and edges, compress, and write GFA."""
    args = parse_args(argv)
    if args.nodes <= 0 or args.length <= 0:
        print(
            "[gen_graph] Node count and length must be positive integers", file=sys.stderr)
        return 2

    rng = random.Random(args.seed)
    sequences = [random_dna(args.length, rng) for _ in range(args.nodes)]
    edges = make_edges(
        num_nodes=args.nodes,
        complexity=args.complexity,
        loop_prob=args.loop_prob,
        rng=rng,
    )
    # Stats before compression
    n0, m0 = len(sequences), len(edges)
    # Perform linear-chain compression
    sequences2, edges2, mapping = compress_linear(sequences, edges)
    n1, m1 = len(sequences2), len(edges2)
    output_path = "-" if args.emit == "stdout" or args.output == "-" else args.output
    write_gfa(output_path, sequences2, edges2)
    # Summary report
    merged_nodes = n0 - n1
    print(
        f"[gen_graph] Original: nodes={n0}, edges={m0}; Compressed: nodes={n1}, edges={m1}; merged_nodes={merged_nodes}")
    print(
        f"[gen_graph] Written: {'stdout' if output_path == '-' else output_path} (len={args.length}, complexity={args.complexity})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
