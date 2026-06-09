# medulla profiling harness

A tiny, dependency-free harness to answer the one question that decides which
optimizations are worth doing:

> **Is medulla CPU-bound (variable/cut evaluation) or I/O-bound (reading +
> decompressing CAF files)?**

The two regimes have completely different "drastic" fixes. Optimizing the wrong
one buys nothing, so **measure before optimizing.**

| If the verdict is... | The lever is... |
|----------------------|-----------------|
| **CPU-bound**        | Confirm the `-O3` build landed (OPTIMIZATIONS.md #1), then the #2 single-pass cut refactor, then micro-opts A–F. |
| **I/O-bound**        | Read fewer bytes: skim/flatten CAFs, prune unused branches, apply a spill/preselection cut, stage files locally instead of streaming xrootd. CPU work (incl. `-O3`) won't move wall time. |

---

## Prerequisites

Run this **inside your medulla runtime environment** (a gpvm / grid node where
ROOT + sbnana are set up) — not on a laptop checkout. Typically:

```bash
source <build>/setup_medulla.sh     # puts the right ROOT/sbnana on the path
```

You need a built `medulla` binary and one real input file referenced by a TOML.
Point at a **single file** (or a short file list) so a run takes seconds–minutes,
not hours.

> ⚠️ Make sure the binary you profile was built optimized. After the CMake
> change (default `Release`/`-O3`), do a **clean** rebuild — a stale build dir
> can silently keep `-O0` and make every measurement meaningless:
> ```bash
> rm -rf build && mkdir build && cd build && cmake .. && make -j
> ```

---

## Quick start

```bash
cd selection/profiling

# 1) Baseline: CPU-vs-I/O from %CPU + cold/warm cache delta (no perf needed).
./profile_medulla.sh --bin /path/to/build/selection/medulla \
                     --toml ../toml/NuMI_nue.toml \
                     --reps 3

# 2) Stronger signal: the branch-count sweep (framework-specific).
python3 make_branch_sweep.py ../toml/NuMI_nue.toml -o /tmp/NuMI_nue_min.toml --keep 2
./profile_medulla.sh --bin /path/to/build/selection/medulla \
                     --toml ../toml/NuMI_nue.toml \
                     --toml-min /tmp/NuMI_nue_min.toml \
                     --reps 3

# 3) Confirmation (Linux + perf available): top hot symbols.
./profile_medulla.sh --bin /path/to/build/selection/medulla \
                     --toml ../toml/NuMI_nue.toml --perf
```

The verdict prints at the end and is saved to `profiling_results/verdict.txt`;
the full log and any perf output land in `profiling_results/`.

---

## The three signals (in order of robustness)

### 1. %CPU  *(needs GNU `/usr/bin/time`; portable)*
The event loop is **single-threaded**. So:
- **%CPU ≈ 100%** → the core is busy computing → **CPU-bound**.
- **%CPU well below 100%** (e.g. 30–60%) with high elapsed → the process is
  waiting on reads → **I/O-bound**.

### 2. Cold vs warm cache  *(portable)*
The first run reads input cold from disk/dCache; repeats hit the OS page cache.
A large **cold/warm wall ratio (>~1.3)** means disk/network read was a big slice
of wall time → I/O-sensitive. (Note: files streamed fresh over xrootd on the
grid are effectively always "cold," so real grid jobs feel the I/O cost even if
a warm local repeat looks fast.)

### 3. Branch-count sweep  *(the framework-specific discriminator)*
Every branch in a tree is built as an independent `SpillMultiVar` that re-runs
the **entire selection cut chain** and rebuilds the particle match-map from
scratch (`framework.cc:708`). A `NuMI_nue.toml` tree books **~113 branches**, so
the selection is recomputed ~113× per interaction.

`make_branch_sweep.py` writes a copy of the *same* selection with only `--keep`
branches per tree. Then:
- **full/reduced wall ≈ 1** → per-branch CPU is noise → **I/O-bound**.
- **full/reduced wall scales up** (≥~1.8×) → per-branch cut re-evaluation
  dominates → **CPU-bound, and OPTIMIZATIONS.md #2 is the lever**.

This is the cleanest test of the #2 hypothesis and needs no perf.

### (optional) perf symbols  *(Linux + perf)*
`--perf` collects `perf stat -dd` and `perf record`/`report`. The script tallies
the heaviest symbols:
- `inflate`/`zlib`/`TBuffer`/`TBranch`/`Unzip`/`xrootd` on top → **I/O/decompress-bound**.
- `pvars::`/`pcuts::`/`vars::`/`std::function` lambdas on top → **CPU-bound**.

---

## Files

| File | Purpose |
|------|---------|
| `profile_medulla.sh`   | Main harness: timing, cold/warm, optional sweep + perf, prints verdict. |
| `make_branch_sweep.py` | Emits a reduced-branch copy of a selection TOML for the sweep. |
| `profiling_results/`   | Logs, `verdict.txt`, perf output (created on first run). |

## Notes / caveats
- Each timed run executes in a fresh temp CWD so the `<analysis>.root` outputs
  don't clobber each other or skew timing.
- Without GNU `time` the script still reports wall time and the cold/warm ratio,
  but cannot compute %CPU — install/locate GNU time or use `--perf`/`--toml-min`
  for a definitive call.
- `make_branch_sweep.py` assumes single-line inline-table branch entries (the
  medulla convention). It warns if it hits a multi-line entry it can't trim.
