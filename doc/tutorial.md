# Tutorial

## Overview

This document provides step-by-step instructions to reproduce the build and the file-conversion and merge experiments reported in the manuscript. It is written to rebuild the project, run the same commands, and reproduce the timing/memory measurements.

**Prerequisites**

- Linux or macOS (experiments were performed on Linux). Windows users may use WSL or a Linux VM.
- Development toolchain: `make`, `gcc`, `g++`, and standard build utilities.
- GNU `time` (the `/usr/bin/time -v` binary) for reproducible resource measurements.

**Repository layout (important paths)**

- Project root: top-level directory containing `Makefile` and `bgfatools` target.
- Source: `src/` (C++ sources and submodules such as `wfa/` and `ksw/`).
- Examples and data used in the paper: `data/graph/` and `data/assembly/`.
- Helper cases used to patch external tools: `case/`.

## Data Preparation

The datasets used in the paper are indicated below:

- 1kcp: [1KCP](https://yanglab.westlake.edu.cn/1kcp/downloadz)
- hprc: [hprc-v1.0-minigraph-chm13.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-chm13.gfa.gz) and [hprc-v1.0-mc-chm13.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13.gfa.gz)
- S.cerevisiae: `./data/assembly/s.c`
- E.coli: [ecoli50.gfa.zst](https://zenodo.org/records/7937947/files/ecoli50.gfa.zst?download=1)

The graphs in the paper were built using two approaches: `minigraph-cactus(mc)` and `minigraph`.

- mc

Repository: https://github.com/ComparativeGenomicsToolkit/cactus

Example command (replace placeholders):

```powershell
cactus-pangenome "${JOBSTORE}" "./${name}.txt" \
                --outDir "${OUT_DIR}" \
                --outName "${name}" \
                --reference "refname" \
                --gfa \
                --maxCores 16 \
                --maxMemory 60G \
                --workDir "${WORK_DIR}"
```

- minigraph

Repository: https://github.com/lh3/minigraph

Example command:

```powershell
minigraph -cxggs file1.fa file2.fa ...
```

Notes: Replace placeholders and ensure fasta inputs are provided in the desired order.

## Build

1. Change into the repository root and build the project:

```powershell
cd bGFA
make clean all
```

Notes about the build:

- The top-level `Makefile` builds the `bgfatools` executable and places intermediate objects in `build/`.
- The default C++ compiler is `g++` (check `CXX` in the `Makefile`). If needed, override with `make CXX=clang++`.
- The `wfa/` subdirectory contains its own Makefile and is built automatically by the top-level `Makefile`.

## Basic commands

All example runs below assume you run them from the repository root and that `bgfatools` is present in `./`.

### Convert

Convert a GFA file to bGFA format (example using S.cerevisiae graph):

```powershell
./bgfatools convert -g ./data/graph/S.C_minigraph.gfa -o ./data/graph/S.C_minigraph.bgfa
```

Measure runtime and peak memory with GNU `time` (example):

```powershell
/usr/bin/time -v ./bgfatools convert -g ./data/graph/S.C_minigraph.gfa -o ./data/graph/S.C_minigraph.bgfa
```

### Merge commands

The repository supports two high-level merge modes. Below are the commands used in the experiments.

Strict merge (type 0):

```powershell
./bgfatools merge -t 0 -g ./file1.gfa ./file2.gfa -o ./merged_strict.bgfa
```

You can also merge already-converted `.bgfa` inputs:

```powershell
./bgfatools merge -t 0 -g ./file1.bgfa ./file2.bgfa -o ./merged_strict.bgfa
```

Exact merge (type 1) — two matching strategies are supported via `-d`:

- Offset-based matching (`-d 1`):

```powershell
./bgfatools merge -t 1 -d 1 -g ./file1.gfa ./file2.gfa -o ./merged_exact_offset.bgfa
```

- Minimizer-based matching (`-d 2`):

```powershell
./bgfatools merge -t 1 -d 2 -g ./file1.gfa ./file2.gfa -o ./merged_exact_minimizer.bgfa
```

Benchmarking and reporting

- Use `/usr/bin/time -v` to capture wall time, user/sys CPU time, and peak resident set size (RSS). For example:

```powershell
/usr/bin/time -v ./bgfatools merge -t 0 -g ./file1.gfa ./file2.gfa -o ./merged_strict.bgfa
```

## Adapting bGFA to external tools (`case/`)

This project includes a `case/` directory with helper source files and patches used to reproduce external tool behavior in the experiments. The following third-party repositories were used and, where noted, were modified by replacing source files with files from `case/`:

- gfatools: https://github.com/lh3/gfatools
- GraphAligner: https://github.com/maickrau/GraphAligner
- Bandage: https://github.com/rrwick/Bandage

If you need to reproduce the modified external tools exactly:

- Extract the files directly
  if a matching `zip` archive is provided under the corresponding directory, extract it directly and build the tool using its original build procedure.
- Or replace the files

  1. Clone the corresponding upstream repository.
  2. Replace the specified source files with the matching files from `case/<toolname>/` (see the `case/` subfolders).
  3. Build the tool following its original build instructions.

Important note about `gfatools`:

Because some replaced files are C++ sources, building the patched `gfatools` version requires using `g++` (or another C++ compiler) instead of `gcc`. When building, ensure the toolchain uses a C++ compiler for C++ sources and links against C++ standard library.

Additional notes

- The commands above were executed and timed on Linux; small differences may appear on other platforms.
- If you encounter build errors, re-run `make clean` and re-run the top-level `make` so the `wfa` and `ksw` objects are rebuilt.
