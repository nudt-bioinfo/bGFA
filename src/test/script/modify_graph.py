'''
file /bGFA/src/test/script/modify_graph.py
brief: Modify at least one node sequence in a GFA
details: Parse S/L/P records from an input GFA, mutate one segment's sequence, and write to a new GFA.
author: grampart
version: 0.1
copyright: Copyright(c) 2025 by grampart, All Rights Reserved.
license: MIT License
'''

import argparse
import os
import random
import sys
from typing import Dict, List, Tuple, Optional

DNA_ALPHABET = "ACGT"


def parse_gfa(path: str) -> Tuple[Dict[int, str], List[Tuple[int, int, str, str, str]], List[Tuple[str, List[str], List[str]]]]:
    """
    Parse a GFA v1 file.
    Returns:
    - segments: map id -> sequence
    - links: list of (from_id, to_id, from_orient, to_orient, overlap)
    - paths: list of (name, node_ids[], offsets[] as strings)
    """
    segments: Dict[int, str] = {}
    links: List[Tuple[int, int, str, str, str]] = []
    paths: List[Tuple[str, List[str], List[str]]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            tag = cols[0]
            if tag == "S":
                # S	<name>	<sequence>
                try:
                    sid = int(cols[1])
                except ValueError:
                    continue
                seq = cols[2] if len(cols) > 2 else ""
                segments[sid] = seq
            elif tag == "L":
                # L	<from>	<from_orient>	<to>	<to_orient>	<overlap>
                if len(cols) < 6:
                    continue
                try:
                    a = int(cols[1])
                    oa = cols[2]
                    b = int(cols[3])
                    ob = cols[4]
                except ValueError:
                    continue
                ov = cols[5]
                links.append((a, b, oa, ob, ov))
            elif tag == "P":
                # P	<path_name>	<segment_names>	<overlaps>
                if len(cols) < 4:
                    continue
                name = cols[1]
                node_ids = cols[2].split(",") if cols[2] else []
                offs = cols[3].split(",") if cols[3] else []
                paths.append((name, node_ids, offs))
    return segments, links, paths


def mutate_sequence(seq: str, rate: float, rng: random.Random) -> str:
    """Mutate sequence with given rate by substituting bases randomly (not equal to original)."""
    if rate <= 0.0:
        return seq
    if rate >= 1.0:
        # Replace entire sequence with random
        return "".join(rng.choice(DNA_ALPHABET) for _ in range(len(seq)))
    out = list(seq)
    k = max(1, int(len(out) * rate))
    idxs = rng.sample(range(len(out)), k)
    for i in idxs:
        orig = out[i]
        choices = [b for b in DNA_ALPHABET if b != orig]
        out[i] = rng.choice(choices)
    return "".join(out)


def write_gfa(path: str, segments: Dict[int, str], links: List[Tuple[int, int, str, str, str]], paths: List[Tuple[str, List[str], List[str]]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("H\tVN:Z:1\n")
        for sid in sorted(segments.keys()):
            f.write(f"S\t{sid}\t{segments[sid]}\n")
        for a, b, oa, ob, ov in links:
            f.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t{ov}\n")
        for name, node_ids, offs in paths:
            nodes_str = ",".join(node_ids) if node_ids else ""
            offs_str = ",".join(offs) if offs else ""
            f.write(f"P\t{name}\t{nodes_str}\t{offs_str}\n")


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Modify node sequences in a GFA: select specific ids or randomly pick N nodes (default 1).")
    p.add_argument("-i", "--input", required=True, help="Input GFA file path")
    p.add_argument("-o", "--output", required=True,
                   help="Output GFA file path")
    p.add_argument("--node", type=int, default=None,
                   help="Single node id to modify (deprecated if --nodes provided)")
    p.add_argument("--nodes", type=str, default=None,
                   help="Comma-separated node ids to modify, e.g. 1,3,5")
    p.add_argument("--count", type=int, default=1,
                   help="Number of random nodes to modify when ids not provided (default 1)")
    p.add_argument("--rate", type=float, default=0.05,
                   help="Mutation rate per node (0..1)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility")
    return p.parse_args(argv)


def _parse_nodes_arg(nodes_arg: Optional[str]) -> Optional[List[int]]:
    if not nodes_arg:
        return None
    try:
        return [int(x) for x in nodes_arg.split(",") if x.strip() != ""]
    except ValueError:
        return None


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    if not os.path.exists(args.input):
        print(f"[modify_graph] Input not found: {args.input}", file=sys.stderr)
        return 2
    rng = random.Random(args.seed)
    segments, links, paths = parse_gfa(args.input)
    if not segments:
        print(
            f"[modify_graph] No segments found in: {args.input}", file=sys.stderr)
        return 3
    # Determine target nodes: from --nodes or --node, else random pick count
    explicit_list = _parse_nodes_arg(args.nodes)
    targets: List[int]
    if explicit_list is not None and len(explicit_list) > 0:
        targets = explicit_list
    elif args.node is not None:
        targets = [args.node]
    else:
        all_ids = sorted(segments.keys())
        k = max(1, min(args.count, len(all_ids)))
        targets = rng.sample(all_ids, k)

    # Validate targets exist
    missing = [t for t in targets if t not in segments]
    if missing:
        print(
            f"[modify_graph] Some target ids not found: {missing}", file=sys.stderr)
        return 4

    # Mutate each target
    changed = 0
    for target in targets:
        orig = segments[target]
        mutated = mutate_sequence(orig, args.rate, rng)
        if mutated == orig:
            mutated = mutate_sequence(
                orig, min(1.0, max(args.rate, 1.0/len(orig) if len(orig) > 0 else 1.0)), rng)
        if mutated != orig:
            changed += 1
        segments[target] = mutated
    write_gfa(args.output, segments, links, paths)
    print(
        f"[modify_graph] Modified nodes: {targets}; changed={changed}; written to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
