"""
Generate GFA graphs with two generation modes.

Mode A (basic):
- Given node count, average length (±20%), and edge multiplier (1.2/1.5/2.0),
  build sequences and a non-linear connection scheme.
- Ensures exactly one source-only node (indeg==0) and one sink-only node (outdeg==0).
  All other nodes have at least one predecessor and one successor.
- The graph is not strictly linear (adds cross edges beyond a simple chain).

Mode B (distribution):
- Given total length and a node-length distribution across 6 bins
  (<10, 10-100, 100-1000, 1k-10k, 10k-100k, >100k),
  either via explicit percentages or presets, estimate node count and sample lengths.
- Then generate connections with the same degree constraints and edge multiplier.

Author: grampart 
Version: 0.1
License: MIT License
"""

from __future__ import annotations

import argparse
import math
import os
import random
import sys
from typing import List, Tuple, Dict, Sequence

DNA_ALPHABET = "ACGT"

# Length bins for mode B: (low, high) inclusive
# For the last bin, the 'high' will be replaced by user-provided max length
BINS_SPECS = [
    (0, 1)
    (1, 9),            # <10
    (10, 100),         # 10-100
    (100, 300),       # 100-1000
    (300, 1000)
    (1000, 3000),     # 1k-10k
    (3000, 10000),     # 1k-10k
    (10000, None),    # >100k (upper bound decided at runtime)
]


PRESET_DISTS = {
    "short": [65, 16, 14, 3, 1, 0.1, 0.02, 0.005],
    "long": [2.5, 7.5, 27.5, 20, 20, 10, 6, 7],
}

EDGE_MULT_CHOICES = {"1.2": 1.2, "1.5": 1.5, "2.0": 2.0}


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA_ALPHABET) for _ in range(max(1, int(length))))


def write_gfa(output_path: str, sequences: Sequence[str],
              edges: Sequence[Tuple[int, int, str, str]]) -> None:
    if output_path == "-":
        f = sys.stdout
        f.write("H\tVN:Z:1\n")
        for i, seq in enumerate(sequences, start=1):
            f.write(f"S\t{i}\t{seq}\n")
        for (a, b, oa, ob) in edges:
            f.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M\n")
        return
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("H\tVN:Z:1\n")
        for i, seq in enumerate(sequences, start=1):
            f.write(f"S\t{i}\t{seq}\n")
        for (a, b, oa, ob) in edges:
            f.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M\n")


def _target_edge_count(n: int, mult: float) -> int:
    # Always at least a chain; then add up to round(mult*n)
    return max(n - 1, int(round(mult * n)))


def _generate_edges(num_nodes: int, edge_mult: float, rng: random.Random) -> List[Tuple[int, int, str, str]]:
    if num_nodes < 2:
        raise ValueError(
            "num_nodes must be >= 2 to ensure unique source and sink")

    # Start with a simple chain to guarantee connectivity and degree constraints for inner nodes
    edges: List[Tuple[int, int, str, str]] = [
        (i, i + 1, "+", "+") for i in range(1, num_nodes)]
    edge_set = set(edges)  # for duplicate checks (includes orientation)

    target = _target_edge_count(num_nodes, edge_mult)

    # Add cross edges forward in index to avoid cycles and preserve source/sink uniqueness
    # Source is 1 (indeg=0), sink is num_nodes (outdeg=0). We avoid edges into 1 or out of N implicitly by a<b.
    while len(edges) < target:
        a = rng.randint(1, num_nodes - 1)
        b = rng.randint(a + 1, num_nodes)
        oa = rng.choice(["+", "-"])
        ob = rng.choice(["+", "-"])
        candidate = (a, b, oa, ob)
        if candidate in edge_set:
            continue
        edge_set.add(candidate)
        edges.append(candidate)

    return edges


def _sample_lengths_mode_a(n: int, avg_len: int, rng: random.Random) -> List[str]:
    if avg_len <= 0:
        raise ValueError("Average length must be positive")
    lo = max(1, int(round(avg_len * 0.8)))
    hi = max(lo, int(round(avg_len * 1.2)))
    seqs = [random_dna(rng.randint(lo, hi), rng) for _ in range(n)]
    return seqs


def _bin_avgs(max_len_last_bin: int) -> List[float]:
    avgs: List[float] = []
    for low, high in BINS_SPECS:
        h = max_len_last_bin if high is None else high
        avgs.append((low + h) / 2.0)
    return avgs


def _normalize_dist(dist: List[float]) -> List[float]:
    s = sum(dist)
    if s <= 0:
        raise ValueError("Distribution sum must be > 0")
    return [x / s for x in dist]


def _perturb_percentages(dist_pct: List[float], jitter: float, rng: random.Random) -> List[float]:
    """Apply ±jitter relative fluctuation to each percentage and rescale to sum 100.

    - jitter=0.1 means each value multiplied by U[0.9, 1.1].
    - Ensures non-negative and final sum exactly 100.0 (within float error).
    """
    jitter = max(0.0, jitter)
    if jitter == 0.0:
        return list(dist_pct)
    perturbed = []
    low = 1.0 - jitter
    high = 1.0 + jitter
    for w in dist_pct:
        factor = rng.uniform(low, high)
        val = max(0.0, w * factor)
        perturbed.append(val)
    total = sum(perturbed)
    if total <= 0:
        return list(dist_pct)
    scale = 100.0 / total
    return [v * scale for v in perturbed]


def _round_counts(total: int, weights: List[float]) -> List[int]:
    # Largest remainder method
    raw = [total * w for w in weights]
    floors = [int(x) for x in raw]
    rem = total - sum(floors)
    fracs = [(i, raw[i] - floors[i]) for i in range(len(weights))]
    fracs.sort(key=lambda t: t[1], reverse=True)
    for i in range(rem):
        floors[fracs[i][0]] += 1
    return floors


def _sample_lengths_mode_b(total_len: int, dist_pct: List[float], max_len_last_bin: int,
                           rng: random.Random, tol: float = 0.1) -> List[str]:
    if total_len <= 0:
        raise ValueError("Total length must be positive")
    if max_len_last_bin < 100001:
        raise ValueError("max length for last bin must be >= 100001")

    weights = _normalize_dist(dist_pct)
    avgs = _bin_avgs(max_len_last_bin)
    exp_avg = sum(w * a for w, a in zip(weights, avgs))
    n_est = max(2, int(round(total_len / max(1.0, exp_avg))))

    counts = _round_counts(n_est, weights)

    # Build ranges with applied max for last bin
    ranges: List[Tuple[int, int]] = []
    for low, high in BINS_SPECS:
        h = max_len_last_bin if high is None else high
        ranges.append((low, h))

    # First pass: sample lengths within their bins; keep track of assigned bin ranges
    seqs: List[str] = []
    lengths: List[int] = []
    assigned_ranges: List[Tuple[int, int]] = []
    for (low, high), cnt in zip(ranges, counts):
        for _ in range(cnt):
            length = rng.randint(low, high)
            lengths.append(length)
            assigned_ranges.append((low, high))
            seqs.append(random_dna(length, rng))

    # In rare cases, if counts rounded down to <2 due to extremely small total_len, pad to 2 nodes
    while len(seqs) < 2:
        seqs.append(random_dna(10, rng))

    # Top-up to reach at least total_len by extending nodes within bin caps, then add nodes if needed
    current_total = sum(lengths)
    # desired bounds for tolerance
    tol = max(0.0, min(1.0, tol))
    lower_bound = int(math.floor(total_len * (1.0 - tol)))
    upper_bound = int(math.ceil(total_len * (1.0 + tol)))

    if current_total < total_len:
        gap = total_len - current_total
        # Extend existing nodes up to their bin upper bounds
        idx = 0
        n_nodes = len(lengths)
        while gap > 0 and idx < n_nodes:
            low, high = assigned_ranges[idx]
            can_inc = max(0, high - lengths[idx])
            if can_inc > 0:
                inc = min(can_inc, gap)
                lengths[idx] += inc
                # append random bases deterministically per node
                seqs[idx] += random_dna(inc, rng)
                gap -= inc
            idx += 1

        # If still not enough, append new nodes sampled by weights
        if gap > 0:
            # Build cumulative weights for sampling bin index
            cum: List[float] = []
            s = 0.0
            for w in weights:
                s += w
                cum.append(s)

            def pick_bin() -> int:
                r = rng.random()
                for i, c in enumerate(cum):
                    if r <= c:
                        return i
                return len(cum) - 1

            while gap > 0:
                bi = pick_bin()
                low, high = ranges[bi]
                # choose a length trying not to overshoot too much, but stay within bin
                # if gap < low, we will overshoot; that's acceptable to satisfy minimum total length
                desired = max(low, min(high, gap))
                seq = random_dna(desired, rng)
                seqs.append(seq)
                lengths.append(desired)
                assigned_ranges.append((low, high))
                # update gap (allow negative to indicate overshoot handled)
                if desired >= gap:
                    gap = 0
                else:
                    gap -= desired

        current_total = sum(lengths)

    # If we exceed the upper bound, try to shrink towards target without going below bin lows
    if current_total > upper_bound:
        # Aim for target total_len if achievable; otherwise get as low as allowed by lows
        target_total = max(lower_bound, min(total_len, current_total))
        shrink = current_total - target_total
        if shrink > 0:
            idx = 0
            n_nodes = len(lengths)
            while shrink > 0 and idx < n_nodes:
                low, _high = assigned_ranges[idx]
                dec_cap = max(0, lengths[idx] - low)
                if dec_cap > 0:
                    dec = min(dec_cap, shrink)
                    lengths[idx] -= dec
                    # trim sequence tail by dec bases
                    if dec > 0:
                        seqs[idx] = seqs[idx][:-
                                              dec] if dec <= len(seqs[idx]) else ""
                    shrink -= dec
                idx += 1

        # Recompute; if still above upper bound, we've hit lower limits; accept best effort
        # If we undershot below lower bound due to rounding, we can extend 1 base back to fix, but unlikely.

    return seqs


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate GFA files with two modes: basic or distribution-driven")

    sub = p.add_subparsers(dest="mode", required=True)

    # Mode A: basic
    pa = sub.add_parser(
        "basic", help="Mode A: specify node count and average length ±20%")
    pa.add_argument("-n", "--nodes", type=int, required=True,
                    help="Number of nodes (>=2)")
    pa.add_argument("-L", "--avg-length", type=int, required=True,
                    help="Average length per node (±20%)")

    # Mode B: distribution
    pb = sub.add_parser(
        "dist", help="Mode B: specify total length and node-length distribution")
    pb.add_argument("--total-length", type=int, required=True,
                    help="Total length across all nodes")
    pb.add_argument("--dist", type=str, default=None,
                    help="Comma-separated percentages for 6 bins (<10,10-100,100-1000,1k-10k,10k-100k,>100k). Sums to 100.")
    pb.add_argument("--preset", type=str, choices=sorted(PRESET_DISTS.keys()), default="long",
                    help="Use a preset node-length distribution")
    pb.add_argument("--max-length", type=int, default=1000000,
                    help=">100k bin upper bound (inclusive), default 1,000,000")
    pb.add_argument("--tolerance", type=float, default=0.10,
                    help="Allowed relative deviation for total length in dist mode (e.g., 0.10 for ±10%)")
    pb.add_argument("--preset-jitter", type=float, default=0.10,
                    help="Relative jitter applied to preset distribution (e.g., 0.10 for ±10% per bin). Set 0 to disable.")

    # Common options
    for x in (pa, pb):
        x.add_argument("--edge-mult", type=str, choices=list(EDGE_MULT_CHOICES.keys()), default="1.5",
                       help="Edge count ≈ edge-mult * nodes (rounded). Choices: 1.2, 1.5, 2.0")
        x.add_argument("-o", "--output", type=str, default="temp/generated.gfa",
                       help="Output GFA path or '-' for stdout")
        x.add_argument("--emit", type=str, choices=["file", "stdout"], default="file",
                       help="Emit target: 'file' or 'stdout' (or use output='-')")
        x.add_argument("--seed", type=int, default=None, help="Random seed")

    return p.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    rng = random.Random(args.seed)
    mult = EDGE_MULT_CHOICES[args.edge_mult]

    if args.mode == "basic":
        if args.nodes < 2:
            print("[gen_gfa] nodes must be >= 2", file=sys.stderr)
            return 2
        sequences = _sample_lengths_mode_a(args.nodes, args.avg_length, rng)
    else:
        # dist mode
        if args.dist is not None:
            try:
                dist = [float(x) for x in args.dist.split(",")]
            except Exception:
                print(
                    "[gen_gfa] --dist must be 6 comma-separated numbers", file=sys.stderr)
                return 2
            if len(dist) != 6:
                print("[gen_gfa] --dist must contain exactly 6 values",
                      file=sys.stderr)
                return 2
        else:
            base = PRESET_DISTS[args.preset]
            # Apply per-bin jitter to presets, then proceed
            dist = _perturb_percentages(
                base, getattr(args, "preset_jitter", 0.10), rng)

            print(dist)
        try:
            sequences = _sample_lengths_mode_b(
                args.total_length, dist, args.max_length, rng, tol=getattr(args, "tolerance", 0.10))
        except ValueError as e:
            print(f"[gen_gfa] {e}", file=sys.stderr)
            return 2

    n = len(sequences)
    try:
        edges = _generate_edges(n, mult, rng)
    except ValueError as e:
        print(f"[gen_gfa] {e}", file=sys.stderr)
        return 2

    output_path = "-" if args.emit == "stdout" or args.output == "-" else args.output
    write_gfa(output_path, sequences, edges)

    # stats
    print(f"[gen_gfa] Nodes={n}, Edges={len(edges)} (mult={mult})")
    if args.mode == "basic":
        print(f"[gen_gfa] Mode=basic, avgLen={args.avg_length}")
    else:
        src = args.dist if args.dist is not None else f"preset:{args.preset}"
        print(
            f"[gen_gfa] Mode=dist, totalLen={args.total_length}, dist={src}, maxLastBin={args.max_length}")
    print(
        f"[gen_gfa] Written: {'stdout' if output_path == '-' else output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
