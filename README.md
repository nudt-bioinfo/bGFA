# bGFA Overview

bGFA is a genome graph toolkit for converting between GFA and BGFA, building minimizer indexes, inspecting graphs, and merging graphs in strict or exact mode.

This document focuses on how to use the commands and how BGFA data is laid out.

Detailed function descriptions and call relationships is in (bGFA-Doxygen)[https://github.com/nudt-bioinfo/bGFA/blob/main/html/index.html]

The commands and execution methods corresponding to the experiments in the paper are described in (bGFA-Tutorial)[https://github.com/nudt-bioinfo/bGFA/blob/main/doc/tutorial.md].

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
