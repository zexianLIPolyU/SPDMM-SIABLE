"""
run_benchmark.py
================

End-to-end state-preparation benchmark: for each ``(n, method)`` cell,
measures both gate counts (single-qubit, CNOT) and wall-clock compile
time, then streams the results into two CSV files used by the plot
scripts in ``plot_figure/``.

Outputs
-------
    data/gate_counts.csv     n, method, u, cx
    data/compile_time.csv    n, method, reps, build_mean, transpile_mean,
                             total_mean, total_ci95   (all times in seconds)

Usage
-----
    python run_benchmark.py
    python run_benchmark.py --nmin 5 --nmax 10
    python run_benchmark.py --methods "PB+BlockZXZ,spdmm,Isometry"
    python run_benchmark.py --reps-overrides "low-rank:2,PB:2"
    python run_benchmark.py --skip-timing      # gate counts only
    python run_benchmark.py --fresh            # wipe CSVs first
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import pandas as pd

from _methods import METHODS, count_gates, time_one, mean_ci


_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, "data")
GATE_CSV = os.path.join(DATA_DIR, "gate_counts.csv")
TIME_CSV = os.path.join(DATA_DIR, "compile_time.csv")


# =============================================================================
# Gate-count experiment
# =============================================================================

def run_gate_counts(n_values, methods, reps, csv_path):
    print("\n[Experiment A: gate counts]")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    first_write = not os.path.exists(csv_path)
    if first_write:
        pd.DataFrame(columns=["n", "method", "u", "cx"]).to_csv(
            csv_path, index=False)

    for n in n_values:
        for method in methods:
            us, cxs = [], []
            for k in range(reps):
                try:
                    u, cx = count_gates(n, k, method)
                except Exception as e:
                    print(f"  n={n} {method} rep={k} FAILED: {e}", flush=True)
                    continue
                us.append(u); cxs.append(cx)
            if not us:
                continue
            u_mean = int(round(sum(us) / len(us)))
            cx_mean = int(round(sum(cxs) / len(cxs)))
            row = {"n": n, "method": method, "u": u_mean, "cx": cx_mean}
            pd.DataFrame([row]).to_csv(csv_path, mode="a", header=False,
                                       index=False)
            print(f"  n={n:2d}  {method:18s}  u={u_mean:>7d}  cx={cx_mean:>7d}",
                  flush=True)


# =============================================================================
# Compile-time experiment
# =============================================================================

def run_compile_time(n_values, methods, reps_default, rep_overrides, csv_path):
    print("\n[Experiment B: compile time]")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    columns = ["n", "method", "reps", "build_mean", "transpile_mean",
               "total_mean", "total_ci95"]
    first_write = not os.path.exists(csv_path)
    if first_write:
        pd.DataFrame(columns=columns).to_csv(csv_path, index=False)

    # Warm-up: compile the smallest case once for each method to amortise
    # any one-time setup cost (Numba JIT, import caches, etc.).
    print("  warming up...")
    for m in methods:
        try:
            time_one(4, 999, m)
        except Exception:
            pass

    for n in n_values:
        for method in methods:
            k = rep_overrides.get(method, reps_default)
            builds, transps, totals = [], [], []
            for r in range(k):
                try:
                    tb, tt, _cx = time_one(n, r, method)
                except Exception as e:
                    print(f"  n={n} {method} rep={r} FAILED: {e}", flush=True)
                    continue
                builds.append(tb); transps.append(tt)
                totals.append(tb + tt)
            if not builds:
                continue
            mb, _ = mean_ci(builds)
            mt, _ = mean_ci(transps)
            mtot, htot = mean_ci(totals)
            row = {"n": n, "method": method, "reps": len(builds),
                   "build_mean": mb, "transpile_mean": mt,
                   "total_mean": mtot, "total_ci95": htot}
            pd.DataFrame([row]).to_csv(csv_path, mode="a", header=False,
                                       index=False)
            print(f"  n={n:2d}  {method:18s}  "
                  f"total = {mtot * 1e3:9.1f} ± {htot * 1e3:6.1f} ms "
                  f"({len(builds)} reps)", flush=True)


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="State-preparation benchmark: gate counts + compile time."
    )
    ap.add_argument("--nmin", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=15)
    ap.add_argument("--reps", type=int, default=5,
                    help="Default reps per compile-time cell.")
    ap.add_argument("--reps-gate", type=int, default=3,
                    help="Reps per gate-count cell.")
    ap.add_argument("--reps-overrides", type=str, default="",
                    help='Per-method compile-time rep overrides, e.g. '
                         '"PB:2,low-rank:2".')
    ap.add_argument("--methods", type=str, default=None,
                    help="Comma-separated subset of methods to benchmark.")
    ap.add_argument("--skip-gates", action="store_true")
    ap.add_argument("--skip-timing", action="store_true")
    ap.add_argument("--fresh", action="store_true",
                    help="Delete existing CSVs before running.")
    ap.add_argument("--append", action="store_true",
                    help="Append to existing CSVs (the default refuses to "
                         "overwrite pre-computed data; use this if you "
                         "actually intended to add rows).")
    args = ap.parse_args()

    n_values = list(range(args.nmin, args.nmax + 1))
    methods = list(METHODS.keys())
    if args.methods:
        methods = [m.strip() for m in args.methods.split(",")]
        for m in methods:
            if m not in METHODS:
                print(f"[error] unknown method: {m}", file=sys.stderr)
                print(f"        available: {list(METHODS.keys())}",
                      file=sys.stderr)
                sys.exit(2)

    rep_overrides = {}
    if args.reps_overrides:
        for kv in args.reps_overrides.split(","):
            k, v = kv.split(":")
            rep_overrides[k.strip()] = int(v)

    if args.fresh:
        for p in (GATE_CSV, TIME_CSV):
            if os.path.exists(p):
                os.remove(p)
                print(f"  removed {p}")
    elif not args.append:
        # Refuse to run if the CSVs already exist, to avoid silently
        # appending to (and thereby polluting) the pre-computed data.
        existing = []
        if not args.skip_gates and os.path.exists(GATE_CSV):
            existing.append(GATE_CSV)
        if not args.skip_timing and os.path.exists(TIME_CSV):
            existing.append(TIME_CSV)
        if existing:
            print("[error] the following CSVs already exist:",
                  file=sys.stderr)
            for p in existing:
                print(f"          {p}", file=sys.stderr)
            print("        Re-run with --fresh to overwrite them, or with "
                  "--append to add rows.", file=sys.stderr)
            sys.exit(3)

    print(f"[config]  n = {n_values[0]}..{n_values[-1]}")
    print(f"[config]  methods: {methods}")
    print(f"[config]  reps_gate = {args.reps_gate}, reps_time = {args.reps}")
    if rep_overrides:
        print(f"[config]  rep overrides: {rep_overrides}")
    print(f"[config]  output CSVs: {GATE_CSV}, {TIME_CSV}")

    t0 = time.perf_counter()
    if not args.skip_gates:
        run_gate_counts(n_values, methods, args.reps_gate, GATE_CSV)
    if not args.skip_timing:
        run_compile_time(n_values, methods, args.reps, rep_overrides, TIME_CSV)
    elapsed = time.perf_counter() - t0

    print(f"\n[done] wall time: {elapsed:.1f} s ({elapsed / 60:.1f} min)")
    print(f"  gate counts  -> {GATE_CSV}")
    print(f"  compile time -> {TIME_CSV}")


if __name__ == "__main__":
    main()
