#!/usr/bin/env python3
"""
make_branch_sweep.py - Produce a reduced-branch copy of a medulla selection TOML
for the CPU-vs-I/O branch-count sweep.

The point of the sweep (see README): run the SAME selection with many branches
vs. only a few. If wall time barely changes, per-branch CPU is noise and the job
is I/O-bound. If wall time scales with branch count, the per-branch cut
re-evaluation (OPTIMIZATIONS.md #2) is the dominant CPU cost.

This does a purely textual trim - it does NOT parse TOML semantics - so it is
robust to includes, custom keys, comments, etc. It relies on the medulla
convention that each branch entry is a single-line inline table inside a
`branch = [ ... ]` array, e.g.:

    branch = [
        {name = "is_nu",          type = "both"},
        {name = "visible_energy", type = "both"},
        ...
    ]

Within every such array it keeps the first --keep entries and drops the rest.
Everything outside branch arrays is copied verbatim.

Usage:
    python3 make_branch_sweep.py NuMI_nue.toml -o NuMI_nue_min.toml --keep 2
"""
import argparse
import sys
import re

# A branch *array opener*: a line whose last non-comment token is `branch = [`.
_OPEN_RE = re.compile(r'^\s*branch\s*=\s*\[\s*(#.*)?$')
# An inline-table entry line, e.g. `{name = "x", type = "both"},`
_ENTRY_RE = re.compile(r'^\s*\{.*\}\s*,?\s*(#.*)?$')
# A closing bracket line (possibly with a trailing comma/comment).
_CLOSE_RE = re.compile(r'^\s*\]\s*,?\s*(#.*)?$')


def reduce_branches(text: str, keep: int) -> tuple[str, int, int]:
    """Return (new_text, n_arrays_trimmed, n_entries_dropped)."""
    out_lines = []
    in_block = False
    kept_in_block = 0
    n_arrays = 0
    n_dropped = 0

    for line in text.splitlines(keepends=True):
        if not in_block:
            out_lines.append(line)
            if _OPEN_RE.match(line):
                in_block = True
                kept_in_block = 0
                n_arrays += 1
            continue

        # Inside a branch array.
        if _CLOSE_RE.match(line):
            out_lines.append(line)
            in_block = False
            continue

        if _ENTRY_RE.match(line):
            if kept_in_block < keep:
                out_lines.append(line)
                kept_in_block += 1
            else:
                n_dropped += 1
            continue

        # Non-entry, non-close line inside the block (blank line, stray
        # comment, or a multi-line entry we don't recognise): keep it verbatim
        # to avoid corrupting the file.
        out_lines.append(line)

    if in_block:
        sys.stderr.write(
            "WARNING: reached EOF while still inside a `branch = [` array; "
            "the input may use multi-line branch entries. Output may be "
            "malformed - inspect it before use.\n"
        )
    return "".join(out_lines), n_arrays, n_dropped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="Path to the full selection TOML.")
    ap.add_argument("-o", "--output", required=True,
                    help="Path to write the reduced-branch TOML.")
    ap.add_argument("--keep", type=int, default=2,
                    help="Branches to keep per tree (default: 2). Keep >=1 so "
                         "the tree is still booked.")
    args = ap.parse_args()

    if args.keep < 1:
        ap.error("--keep must be >= 1 (a tree with zero branches books nothing).")

    try:
        with open(args.input, "r") as fh:
            text = fh.read()
    except OSError as e:
        sys.stderr.write(f"ERROR: cannot read {args.input}: {e}\n")
        return 1

    new_text, n_arrays, n_dropped = reduce_branches(text, args.keep)

    if n_arrays == 0:
        sys.stderr.write(
            "ERROR: found no `branch = [` arrays in the input. Is this a "
            "selection TOML? (If branches are pulled in via `include`, run "
            "this on the included file instead.)\n"
        )
        return 1

    try:
        with open(args.output, "w") as fh:
            fh.write(new_text)
    except OSError as e:
        sys.stderr.write(f"ERROR: cannot write {args.output}: {e}\n")
        return 1

    print(f"Wrote {args.output}")
    print(f"  trimmed {n_arrays} branch array(s), keeping {args.keep} each")
    print(f"  dropped {n_dropped} branch entries total")
    print("\nNow run, e.g.:")
    print(f"  ./profile_medulla.sh --bin <medulla> --toml {args.input} "
          f"--toml-min {args.output} --reps 3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
