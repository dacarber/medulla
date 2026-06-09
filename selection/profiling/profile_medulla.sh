#!/usr/bin/env bash
#
# profile_medulla.sh - Determine whether medulla is CPU-bound or I/O-bound.
#
# This harness answers the one question that decides which optimizations are
# worth doing:  is wall time dominated by CPU (variable/cut evaluation) or by
# I/O (reading + decompressing CAF files)?  The two regimes have completely
# different "drastic" fixes, so measure before optimizing.
#
# It uses three signals, in order of robustness (no perf / no recompile needed
# for the first two):
#
#   1. %CPU            - single-threaded loop pinned near 100% => CPU-bound;
#                        well below 100% with high elapsed     => I/O/wait-bound.
#   2. cold vs warm    - first run (cold page cache) vs repeats (warm). A large
#                        cold>warm gap means the OS page cache / disk read was a
#                        big fraction of wall time => I/O-sensitive.
#   3. perf (optional) - top symbols: ROOT/zlib/TBuffer/inflate dominating
#                        => I/O/decompress-bound;  vars::/pcuts::/pvars:: lambdas
#                        dominating => CPU-bound.
#
# Optional branch-count sweep (--toml-min): time a full-branch TOML against a
# few-branch copy of the SAME selection. If wall time barely moves, the per-
# branch CPU is noise (I/O-bound). If it scales with branch count, the
# per-branch cut re-evaluation (OPTIMIZATIONS.md #2) is the lever.
#
# Run this INSIDE your medulla environment:
#     source <build>/setup_medulla.sh
#     ./profile_medulla.sh --bin <build>/selection/medulla --toml <cfg.toml>
#
# Author: profiling harness for medulla
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults / argument parsing
# ---------------------------------------------------------------------------
BIN=""
TOML=""
TOML_MIN=""
REPS=3
USE_PERF=0
OUTDIR="profiling_results"

usage() {
    cat <<'EOF'
Usage: profile_medulla.sh --bin <medulla> --toml <cfg.toml> [options]

Required:
  --bin   PATH    Path to the compiled `medulla` executable.
  --toml  PATH    Selection TOML to run (the full, real workload).

Options:
  --reps  N       Number of timed repetitions after warm-up (default: 3).
  --toml-min PATH A copy of the same selection with only a few branches, for
                  the branch-count sweep (see README). Enables the CPU-vs-I/O
                  slope test.
  --perf          Also collect `perf stat` + `perf record/report` (Linux only,
                  requires perf to be installed and usable).
  --out   DIR     Output directory for logs (default: profiling_results).
  -h, --help      Show this help.

The verdict is printed at the end and written to <out>/verdict.txt.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bin)      BIN="$2"; shift 2 ;;
        --toml)     TOML="$2"; shift 2 ;;
        --toml-min) TOML_MIN="$2"; shift 2 ;;
        --reps)     REPS="$2"; shift 2 ;;
        --perf)     USE_PERF=1; shift ;;
        --out)      OUTDIR="$2"; shift 2 ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$BIN" || -z "$TOML" ]]; then
    echo "ERROR: --bin and --toml are required." >&2
    usage; exit 1
fi
if [[ ! -x "$BIN" ]]; then
    echo "ERROR: medulla binary not found or not executable: $BIN" >&2
    exit 1
fi
if [[ ! -f "$TOML" ]]; then
    echo "ERROR: TOML not found: $TOML" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$OUTDIR/profile_${STAMP}.log"
VERDICT="$OUTDIR/verdict.txt"

# ---------------------------------------------------------------------------
# Environment detection
# ---------------------------------------------------------------------------
GNU_TIME=""
if command -v /usr/bin/time >/dev/null 2>&1 && /usr/bin/time -f "%e" true >/dev/null 2>&1; then
    GNU_TIME="/usr/bin/time"     # GNU time supports -f / -v
fi
NPROC="$( (command -v nproc >/dev/null 2>&1 && nproc) || getconf _NPROCESSORS_ONLN 2>/dev/null || echo "?")"
HAVE_PERF=0
if [[ "$USE_PERF" -eq 1 ]]; then
    if command -v perf >/dev/null 2>&1 && perf stat true >/dev/null 2>&1; then
        HAVE_PERF=1
    else
        echo "WARNING: --perf requested but perf is unavailable/unusable; skipping perf." | tee -a "$LOG"
    fi
fi

log() { echo "$@" | tee -a "$LOG"; }

log "============================================================"
log " medulla profiling harness   ($STAMP)"
log "============================================================"
log "binary    : $BIN"
log "toml      : $TOML"
log "reps      : $REPS"
log "cores     : $NPROC (note: the event loop is single-threaded)"
log "GNU time  : ${GNU_TIME:-<not found - using shell builtin fallback>}"
log "perf      : $([[ $HAVE_PERF -eq 1 ]] && echo yes || echo no)"
log ""

# ---------------------------------------------------------------------------
# Timing helpers
#   run_timed <label> <toml>  ->  echoes "wall_s user_s sys_s pct_cpu maxrss_kb"
# ---------------------------------------------------------------------------
run_timed() {
    local label="$1" toml="$2"
    local tf; tf="$(mktemp)"
    # medulla writes <analysis-name>.root in CWD; run in a scratch dir so reps
    # don't clobber each other and so output writes don't skew the measurement.
    local scratch; scratch="$(mktemp -d)"
    if [[ -n "$GNU_TIME" ]]; then
        # %e wall, %U user, %S sys, %P percent-cpu, %M maxRSS(kB)
        ( cd "$scratch" && "$GNU_TIME" -f "%e %U %S %P %M" -o "$tf" \
            "$BIN" "$toml" >/dev/null 2>>"$LOG" ) || { echo "RUN FAILED ($label) - see $LOG" >&2; rm -rf "$scratch"; return 1; }
        # %P looks like "97%" - strip the percent sign.
        awk '{gsub(/%/,"",$4); print $1, $2, $3, $4, $5}' "$tf"
    else
        # Portable fallback: wall via date, user+sys via bash `times`, no RSS.
        local t0 t1 wall
        t0="$(date +%s.%N)"
        ( cd "$scratch" && "$BIN" "$toml" >/dev/null 2>>"$LOG" ) || { echo "RUN FAILED ($label) - see $LOG" >&2; rm -rf "$scratch"; return 1; }
        t1="$(date +%s.%N)"
        wall="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')"
        # No reliable per-child CPU split without GNU time; report wall only.
        echo "$wall NA NA NA NA"
    fi
    rm -f "$tf"; rm -rf "$scratch"
}

# ---------------------------------------------------------------------------
# 1) Warm-up (prime page cache + dynamic linker), then timed reps
# ---------------------------------------------------------------------------
log "------------------------------------------------------------"
log " Step 1: warm-up + $REPS timed reps on the full TOML"
log "------------------------------------------------------------"
log "Cold run (page cache cold; first touch of input files)..."
COLD="$(run_timed cold "$TOML")"
COLD_WALL="$(echo "$COLD" | awk '{print $1}')"
log "  cold:  wall=${COLD_WALL}s   (user/sys/%cpu/maxRSS: $(echo "$COLD" | awk '{print $2"/"$3"/"$4"/"$5}'))"

WARM_WALLS=()
WARM_PCTS=()
WARM_USER=()
WARM_SYS=()
MAXRSS=0
for ((r=1; r<=REPS; r++)); do
    OUT="$(run_timed "warm$r" "$TOML")"
    w="$(echo "$OUT" | awk '{print $1}')"
    u="$(echo "$OUT" | awk '{print $2}')"
    s="$(echo "$OUT" | awk '{print $3}')"
    p="$(echo "$OUT" | awk '{print $4}')"
    m="$(echo "$OUT" | awk '{print $5}')"
    WARM_WALLS+=("$w"); WARM_PCTS+=("$p"); WARM_USER+=("$u"); WARM_SYS+=("$s")
    [[ "$m" =~ ^[0-9]+$ ]] && (( m > MAXRSS )) && MAXRSS=$m
    log "  warm$r: wall=${w}s   user=${u}s sys=${s}s  %cpu=${p}  maxRSS=${m}kB"
done

# Aggregate warm stats (median wall, mean %cpu where available).
WARM_WALL_MED="$(printf '%s\n' "${WARM_WALLS[@]}" | sort -n | awk '{a[NR]=$1} END{print (NR%2)? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2}')"
PCT_MEAN="$(printf '%s\n' "${WARM_PCTS[@]}" | awk '($1!="NA"){s+=$1;n++} END{if(n) printf "%.0f", s/n; else print "NA"}')"
USER_MEAN="$(printf '%s\n' "${WARM_USER[@]}" | awk '($1!="NA"){s+=$1;n++} END{if(n) printf "%.2f", s/n; else print "NA"}')"
SYS_MEAN="$(printf '%s\n' "${WARM_SYS[@]}" | awk '($1!="NA"){s+=$1;n++} END{if(n) printf "%.2f", s/n; else print "NA"}')"

# Cold/warm ratio (I/O sensitivity).
COLD_WARM_RATIO="$(awk -v c="$COLD_WALL" -v w="$WARM_WALL_MED" 'BEGIN{ if(w>0) printf "%.2f", c/w; else print "NA"}')"

log ""
log "  median warm wall : ${WARM_WALL_MED}s"
log "  mean %CPU        : ${PCT_MEAN}"
log "  mean user / sys  : ${USER_MEAN}s / ${SYS_MEAN}s"
log "  cold/warm ratio  : ${COLD_WARM_RATIO}  (>1.3 suggests I/O-sensitive)"
log "  peak RSS         : ${MAXRSS}kB"
log ""

# ---------------------------------------------------------------------------
# 2) Optional branch-count sweep (the framework-specific discriminator)
# ---------------------------------------------------------------------------
SWEEP_VERDICT=""
if [[ -n "$TOML_MIN" ]]; then
    if [[ ! -f "$TOML_MIN" ]]; then
        log "WARNING: --toml-min '$TOML_MIN' not found; skipping sweep."
    else
        log "------------------------------------------------------------"
        log " Step 2: branch-count sweep (full vs reduced-branch TOML)"
        log "------------------------------------------------------------"
        run_timed warmup-min "$TOML_MIN" >/dev/null 2>&1 || true
        MIN_OUT="$(run_timed min "$TOML_MIN")"
        MIN_WALL="$(echo "$MIN_OUT" | awk '{print $1}')"
        log "  reduced-branch wall : ${MIN_WALL}s"
        log "  full-branch wall    : ${WARM_WALL_MED}s"
        SWEEP_DELTA="$(awk -v f="$WARM_WALL_MED" -v m="$MIN_WALL" 'BEGIN{ if(m>0) printf "%.2f", f/m; else print "NA"}')"
        log "  full/reduced ratio  : ${SWEEP_DELTA}"
        # If the full (many-branch) run is much slower than the reduced one,
        # per-branch CPU dominates => CPU-bound and #2 is the lever.
        SWEEP_VERDICT="$(awk -v d="$SWEEP_DELTA" 'BEGIN{
            if(d=="NA"){print "inconclusive"}
            else if(d>=1.8){print "CPU-bound (per-branch evaluation dominates -> OPTIMIZATIONS.md #2)"}
            else if(d<=1.2){print "I/O-bound (branch count barely changes wall time)"}
            else {print "mixed"}
        }')"
        log "  sweep verdict       : ${SWEEP_VERDICT}"
        log ""
    fi
fi

# ---------------------------------------------------------------------------
# 3) Optional perf
# ---------------------------------------------------------------------------
PERF_VERDICT=""
if [[ "$HAVE_PERF" -eq 1 ]]; then
    log "------------------------------------------------------------"
    log " Step 3: perf stat + perf record/report"
    log "------------------------------------------------------------"
    PERF_STAT="$OUTDIR/perf_stat_${STAMP}.txt"
    PERF_DATA="$OUTDIR/perf_${STAMP}.data"
    PERF_REPORT="$OUTDIR/perf_report_${STAMP}.txt"
    scratch="$(mktemp -d)"
    ( cd "$scratch" && perf stat -dd -o "$PERF_STAT" "$BIN" "$TOML" >/dev/null 2>>"$LOG" ) || log "  perf stat failed (see $LOG)"
    ( cd "$scratch" && perf record -F 999 -g -o "$PERF_DATA" "$BIN" "$TOML" >/dev/null 2>>"$LOG" ) || log "  perf record failed (see $LOG)"
    rm -rf "$scratch"
    if [[ -f "$PERF_DATA" ]]; then
        perf report -i "$PERF_DATA" --stdio 2>>"$LOG" | head -40 > "$PERF_REPORT" || true
        log "  perf stat   -> $PERF_STAT"
        log "  perf report -> $PERF_REPORT (top 40 lines)"
        # Classify the heaviest symbols.
        IO_HITS="$(grep -ciE 'inflate|zlib|libz|TBuffer|TBranch|Unzip|R__unzip|decompress|GetEntry|TTree|xrootd|XrdCl' "$PERF_REPORT" || true)"
        CPU_HITS="$(grep -ciE 'pvars|pcuts|vars::|cuts::|spill_multivar|std::function|SRProxy|operator\(\)' "$PERF_REPORT" || true)"
        log "  top-symbol tally: I/O-ish=${IO_HITS}  CPU-ish=${CPU_HITS}"
        PERF_VERDICT="$(awk -v io="$IO_HITS" -v cpu="$CPU_HITS" 'BEGIN{
            if(io>cpu){print "I/O/decompress-bound (ROOT/zlib symbols on top)"}
            else if(cpu>io){print "CPU-bound (medulla var/cut lambdas on top)"}
            else {print "mixed/inconclusive"}
        }')"
        log "  perf verdict     : ${PERF_VERDICT}"
    fi
    log ""
fi

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
{
    echo "============================================================"
    echo " VERDICT  ($STAMP)"
    echo "============================================================"
    echo "median warm wall : ${WARM_WALL_MED}s"
    echo "mean %CPU        : ${PCT_MEAN} (single-thread; ~100% => CPU-bound)"
    echo "cold/warm ratio  : ${COLD_WARM_RATIO} (>1.3 => I/O-sensitive)"
    [[ -n "$SWEEP_VERDICT" ]] && echo "branch sweep     : ${SWEEP_VERDICT}"
    [[ -n "$PERF_VERDICT"  ]] && echo "perf symbols     : ${PERF_VERDICT}"
    echo ""

    # Primary heuristic from %CPU (most robust portable signal).
    primary="inconclusive"
    if [[ "$PCT_MEAN" != "NA" ]]; then
        if   awk -v p="$PCT_MEAN" 'BEGIN{exit !(p>=85)}'; then primary="CPU-bound"
        elif awk -v p="$PCT_MEAN" 'BEGIN{exit !(p<=60)}'; then primary="I/O-bound"
        else primary="mixed"; fi
    fi
    echo "PRIMARY SIGNAL (%CPU): ${primary}"
    echo ""
    echo "Recommended next step:"
    case "$primary" in
        CPU-bound)
            echo "  * CPU dominates. First CONFIRM the -O3 build actually took effect"
            echo "    (clean rebuild; OPTIMIZATIONS.md #1). Then the big code lever is"
            echo "    #2: evaluate the selection cut + particle match-map ONCE per event"
            echo "    instead of once per branch (~Nx fewer cut evaluations)."
            echo "  * Micro-opts A-F in OPTIMIZATIONS.md compound on top."
            ;;
        I/O-bound)
            echo "  * I/O dominates. CPU optimizations (incl. -O3) will NOT move wall time"
            echo "    much. Focus on reading fewer bytes: skim/flatten CAFs, prune unused"
            echo "    branches, apply a spill/preselection cut, and stage files locally"
            echo "    instead of streaming over xrootd (OPTIMIZATIONS.md #3-style I/O work)."
            ;;
        mixed)
            echo "  * Mixed regime. Do BOTH: confirm -O3, then the #2 single-pass cut"
            echo "    refactor for CPU, plus branch-pruning/skimming for I/O."
            ;;
        *)
            echo "  * Could not classify from %CPU (no GNU time?). Re-run with GNU"
            echo "    /usr/bin/time available, or use --perf, or the --toml-min sweep."
            ;;
    esac
    echo "============================================================"
} | tee "$VERDICT" | tee -a "$LOG"

echo ""
echo "Full log : $LOG"
echo "Verdict  : $VERDICT"
