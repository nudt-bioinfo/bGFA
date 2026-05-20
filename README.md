# bGFA Overview

bGFA is a genome graph toolkit for converting between GFA and BGFA, building minimizer indexes, inspecting graphs, and merging graphs in strict or exact mode.

This document focuses on how to use the commands and how BGFA data is laid out.

## Fast Start

### Compilation

`make clean all`

### Usage

- `convert`: Convert `GFA <-> BGFA`.
- `index`: Build a `.bmin` minimizer index from `GFA` or `BGFA`.
- `view`: Print graph or minimizer information.
- `merge`: Merge one or more graphs into a single output.
- `stats`: Show graph statistics.

## Command Guide

### `convert`

Use this when you want to switch file formats.

```bash
./bgfatools convert -g input.gfa -o auto
./bgfatools convert -g input.bgfa -o auto
```

Common options:

- `-g, --gfa`: Input file. The command accepts both `.gfa` and `.bgfa`.
- `-o, --output`: Output file. `auto` uses the opposite suffix automatically.
- `--output-mode all|split`: Choose whether to write a single BGFA file or split output.
- `--path-mode direct|inc`: Choose path encoding mode.
- `--normalized`: Normalize linear one-in-one-out segments once after loading GFA.
- `--segment-no-id`: Store BGFA segments without explicit segment ID.

When to choose which mode:

- Use `all` for normal conversion.
- Use `split` when you want split BGFA outputs for downstream processing.
- Use `inc` when you prefer compact path encoding.
- Use `direct` when you want direct node encoding in paths.

### `index`

Use this to build a minimizer index for later inspection or graph matching.

```bash
./bgfatools index -g input.gfa -o auto
./bgfatools index -b input.bgfa -o auto
```

Common options:

- `-g, --gfa` / `-b, --bgfa`: Pick the input format.
- `-m, --mode`: Currently use `minimizer` for indexing.
- `--minimizer-mode single|single-reverse|double`: Choose how strand handling works.
- `-k, --kmer-length`: Minimizer k-mer length.
- `-w, --window`: Minimizer window size.
- `-o, --output`: Output `.bmin` path.

Selection hints:

- Use `double` for the most complete indexing behavior.
- Use `single-reverse` if you want a lighter mode that still considers reverse complements.
- Use `single` only when you want forward-strand traversal only.

### `view`

Use this to inspect a graph or minimizer index.

```bash
./bgfatools view -g merged.bgfa
./bgfatools view -m sample.bmin
```

- `-g, --bgfa`: Print graph content.
- `-m, --bmin`: Print minimizer index content.

### `merge`

Use this to combine multiple graphs.

```bash
./bgfatools merge -g a.bgfa -g b.bgfa -o merged.bgfa -t 0
./bgfatools merge -g a.bgfa -g b.bgfa -o merged.bgfa -t 1 -l 90 -k wfa
```

Common options:

- `-g, --graph`: One or more input graph files.
- `-o, --output`: Output merged file.
- `-t, --type`: `0` for strict merge, `1` for exact merge.
- `-l, --merge-limited`: Exact-merge similarity threshold.
- `-k, --alignment-kernel`: Alignment backend, `ksw` or `wfa`.
- `-d, --details-type`: Extra exact-merge aid: `0` none, `1` distance, `2` minimizer.

Selection hints:

- Use `strict` when graphs should only merge on direct consistency.
- Use `exact` when sequence-level alignment and splitting are needed.
- Use `distance` if graph distance can quickly filter candidates.
- Use `minimizer` if you want minimizer screening before alignment.
- Use `wfa` for the default exact-merge backend; use `ksw` when you want the alternative backend.

### `stats`

Use this to print summary statistics for a graph.

```bash
./bgfatools stats -g merged.bgfa
```

## Directory Snapshot

- `src/`: Main implementation for graph storage, conversion, indexing, merge, alignment, and CLI.
- `test/`: Regression and behavior tests.
- `case/`: Adapted code for Bandage, gfatools, and GraphAligner.
- `build/`: Optional build artifacts.
- `data/`: Example inputs and outputs.

## `case/` Adaptation Notes

The `case/` directory contains code adapted from external projects so they can work with the current bGFA library and tools. These files are intended to replace the original project files after adaptation.

- `case/Bandage/`: Bandage-side adaptation files such as `assemblygraph.cpp`, `assemblygraph.h`, `bgfa.hpp`, and `bgfa_graph.hpp`.
- `case/gfatools/`: gfatools-side adaptation files such as `gfa-io.c`, `bgfa.hpp`, and `bgfa_graph.hpp`.
- `case/GraphAligner/`: GraphAligner-side adaptation files such as `Aligner.cpp`, `AlignerMain.cpp`, `GfaGraph.cpp`, `GfaGraph.h`, `bgfa.hpp`, and `bgfa_graph.hpp`.

## `src/` Layout

The main source directory is organized around the CLI and graph pipeline:

- `main.cpp`: CLI entry point.
- `bgfa_args.hpp`: Argument structures and option mappings.
- `bgfa_subcommand.hpp`: Command handlers.
- `bgfa_graph.hpp`: Graph, segment, link, path, and walk data structures.
- `bgfa_merge.hpp`: Merge logic.
- `bgfa_alignment.hpp`: Alignment wrappers.
- `bgfa_seed.hpp`: Minimizer indexing and `.bmin` IO.
- `compress.hpp`: Packing helpers.
- `ksw/`, `wfa/`: Vendored alignment dependencies.

## BGFA Binary Layout

The BGFA file starts with a small header, then stores segment data, link data, path data, and optional walk data.

### Header

The first words record offsets so the reader can jump to each block quickly.

```text
pre_info    -> flags for path mode and segment-id storage
s_offset    -> start of segment block
l_offset    -> start of link block
p_offset    -> start of path block
w_offset    -> start of walk block
```

### Segment block

Segment records store one graph node's sequence and a small amount of metadata. In practice, a segment is the binary form of `Segment` in `src/bgfa_graph.hpp`.

What is stored for each segment:

- `id`: segment ID.
- `dis`: graph distance used by some merge and screening paths.
- `length`: sequence length.
- `sequence`: bases encoded in binary form.
- optional storage flags: whether the segment ID is written explicitly and how the sequence is packed.

How to understand the packed content:

- The sequence itself is encoded with the internal base mapping `A/C/G/T -> 0/1/2/3`.
- Short segments may use a compact representation to save space.
- Longer segments use a full representation with explicit length information and packed base words.
- The reader and writer follow the same segment packing contract through `GFA::write2Bgfa` and `GFA::loadFromFileBgfa`.

Practical reading tips:

- Treat the segment block as the canonical node table of the graph.
- If you only need graph topology, segment IDs and lengths are usually enough to inspect first.
- If you need sequence validation, read the stored binary bases back to a string and compare the original sequence.
- When `--segment-no-id` is enabled in `convert`, the segment ID is not stored explicitly in the record and is reconstructed by read order.

In short, the segment block is the place where BGFA keeps the graph's node sequences, lengths, and the metadata needed to rebuild the graph later.

### Link block

Links are stored as a compressed adjacency list.

```text
link_num
n_rows
indptr[0..n_rows]
link_values...
```

Interpretation:

- `link_num`: total number of links.
- `n_rows`: number of source nodes with outgoing links.
- `indptr`: prefix sum array that marks the start of each source node's link list.
- `link_values`: each link is packed as `(to_id << 2) | dir_bits`.

Direction bits:

- Low 2 bits encode the edge orientation.
- `dir_bits = (from_dir << 1) | to_dir`.
- `dir` values are mapped from `+` and `-` through `GlobalVariant::dir2bin`.

Practical reading rule:

- Shift right by `2` to get the target segment ID.
- Mask with `0x03` to get the direction bits.

### Path block

`Path` stores an ordered node list.

```text
n_nodes
node_0
node_1
...
node_(n-1)
```

Notes:

- Each node is encoded as a `uword_t` value.
- In direct mode, a node is stored as `(node_id << 1) | dir_bit`.
- In increment mode, the internal node list can use a compact increment-based representation.
- When converting to text, the path is printed as an ordered sequence of oriented nodes.

### Walk block

`Walk` extends `Path` with sample-level metadata.

```text
sample_id
haplotype_id
sequence_id
path data...
```

Notes:

- Walk records are useful when the graph keeps sample, haplotype, or sequence identifiers.
- The embedded path is serialized immediately after the three IDs.
- Walk text output is formatted like a GFA `W` record.

### Merging-related flags

These high-bit flags are used internally when graphs are merged:

- `MERGE_PLACEHOLDER = 4`: used for graph ID shifting.
- `EXACT_MERGE_PLACEHOLDER = 12`: used for exact-merge renumbering.
- `EXACT_MERGE_CLEAR_MASK`: clears exact-merge tags from stored node values.

## Binary Reading/Writing

The main entry points are:

- `GFA::write2Bgfa`: Write the BGFA binary file.
- `GFA::loadFromFileBgfa`: Read the BGFA binary file.
- `GFA::write2GFA`: Export a text GFA file.
- `GFA::write2BgfaSplit`: Write split BGFA output when needed.

## Picking the Right Command

- Use `convert` first if your input is not already in the format you want.
- Use `index` when you need a `.bmin` for inspection or follow-up tools.
- Use `merge -t 0` for conservative graph combining.
- Use `merge -t 1` when sequence-level alignment is required.
- Use `view` before and after conversion to confirm the graph still looks correct.
- Use `stats` when you only need a quick summary.

## Notes

- Build with `bgfatools` or rename the output as needed for your workflow.
- Merge behavior depends on the chosen type, threshold, and backend.
- The `case/` folder is for adapting bGFA to other projects, not for core library changes.
- Similarity score tests live in `test/test_SimilarScore.cpp` and `test/README.md`.

(This README now emphasizes usage and file layout rather than implementation details.)
