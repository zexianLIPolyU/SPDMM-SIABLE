"""
run_siable_benchmark.py
=======================

End-to-end benchmark for the SIABLE block-encoding algorithm.  Three
experiments are performed in a single pass:

  (A) Full-rank SIABLE -- CNOT and single-qubit gate counts vs n
  (B) Low-rank SIABLE  -- CNOT counts at fixed n, varying rank K,
                          with the auto-dispatch method recorded
  (C) Full-rank SIABLE -- compile time (build + transpile) vs n

The matrix sizes used are ``2^(n-1) x 2^(n-1)``, so the total circuit
has ``n`` qubits including the ancilla (paper convention).

Outputs (under ``data/``):
    siable_full_rank.csv     n, u, cx, build_mean, transpile_mean,
                             total_mean, total_ci95
    siable_low_rank.csv      n, K, method, u, cx, predicted_full_rank,
                             predicted_isometry, predicted_spdmm_half,
                             predicted_state_prep

Usage
-----
    python run_siable_benchmark.py
    python run_siable_benchmark.py --nmax 6 --reps 3
    python run_siable_benchmark.py --skip-low-rank
    python run_siable_benchmark.py --skip-full-rank
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from qiskit import transpile

# Make the local 'spdmm' package importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from spdmm import siable, siable_low_rank  # noqa: E402

# scipy.stats only used for the Student-t CI helper
from _methods import mean_ci  # noqa: E402


DATA_DIR = os.path.join(_HERE, "data")
FULL_CSV = os.path.join(DATA_DIR, "siable_full_rank.csv")
LOW_CSV  = os.path.join(DATA_DIR, "siable_low_rank.csv")


# =============================================================================
# Random matrix generators
# =============================================================================

def _random_matrix(N, seed, scale=0.3):
    """Random complex matrix of size N x N with spectral norm <= 1."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    A *= scale
    # Be defensive about spectral norm
    s = np.linalg.svd(A, compute_uv=False)
    if s[0] > 1.0:
        A /= s[0]
    return A


def _random_rank_K_matrix(N, K, seed):
    """Random N x N matrix of exact rank K (W @ V with W: N x K, V: K x N)."""
    rng = np.random.default_rng(seed)
    W = rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))
    V = rng.standard_normal((K, N)) + 1j * rng.standard_normal((K, N))
    A = W @ V
    nrm = max(np.linalg.norm(A, 2), 1.0)
    return A / nrm


# =============================================================================
# (A) Full-rank SIABLE: gate counts + compile time vs n
# =============================================================================

def run_full_rank(n_values, reps, csv_path):
    print("\n[Experiment A+C: full-rank SIABLE, gate counts + compile time]")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    columns = ["n", "u", "cx", "reps", "build_mean", "transpile_mean",
               "total_mean", "total_ci95"]
    if not os.path.exists(csv_path):
        pd.DataFrame(columns=columns).to_csv(csv_path, index=False)

    for n in n_values:
        N = 2 ** (n - 1)         # data dimension; ancilla is the (n-th) qubit
        A = _random_matrix(N, seed=42 + n)
        builds, transps, totals = [], [], []
        u_count, cx_count = None, None
        for r in range(reps):
            try:
                t0 = time.perf_counter()
                qc, alpha, info = siable(A, return_info=True)
                t_build = time.perf_counter() - t0
                t0 = time.perf_counter()
                tc = transpile(qc, basis_gates=["u", "cx"],
                               optimization_level=0)
                t_transpile = time.perf_counter() - t0
            except Exception as e:
                print(f"  n={n} rep={r} FAILED: {e}", flush=True)
                continue
            builds.append(t_build)
            transps.append(t_transpile)
            totals.append(t_build + t_transpile)
            ops = tc.count_ops()
            u_count = int(ops.get("u", 0))
            cx_count = int(ops.get("cx", 0))
        if not builds:
            continue
        mb, _ = mean_ci(builds)
        mt, _ = mean_ci(transps)
        mtot, htot = mean_ci(totals)
        row = {"n": n, "u": u_count, "cx": cx_count, "reps": len(builds),
               "build_mean": mb, "transpile_mean": mt,
               "total_mean": mtot, "total_ci95": htot}
        pd.DataFrame([row]).to_csv(csv_path, mode="a", header=False,
                                   index=False)
        print(f"  n={n:2d}  u={u_count:>7d}  cx={cx_count:>7d}  "
              f"total = {mtot * 1e3:9.1f} ± {htot * 1e3:6.1f} ms "
              f"({len(builds)} reps)", flush=True)


# =============================================================================
# (B) Low-rank SIABLE: which method does auto-dispatch select?
# =============================================================================

def run_low_rank(n_values, K_values, csv_path):
    print("\n[Experiment B: low-rank SIABLE, method dispatch]")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    columns = ["n", "K", "method", "u", "cx",
               "predicted_state_prep", "predicted_isometry",
               "predicted_spdmm_half", "predicted_full_rank"]
    if not os.path.exists(csv_path):
        pd.DataFrame(columns=columns).to_csv(csv_path, index=False)

    for n in n_values:
        N = 2 ** (n - 1)
        for K in K_values:
            if K > N:
                continue
            try:
                A = _random_rank_K_matrix(N, K, seed=100 + n * 17 + K)
                qc, alpha, info = siable_low_rank(A, K, return_info=True)
                tc = transpile(qc, basis_gates=["u", "cx"],
                               optimization_level=0)
            except Exception as e:
                print(f"  n={n} K={K} FAILED: {e}", flush=True)
                continue
            ops = tc.count_ops()
            u_count = int(ops.get("u", 0))
            cx_count = int(ops.get("cx", 0))
            pred = info.get("predicted_cost", {})
            row = {"n": n, "K": K, "method": info["method"],
                   "u": u_count, "cx": cx_count,
                   "predicted_state_prep": pred.get("state_prep", ""),
                   "predicted_isometry":   pred.get("isometry", ""),
                   "predicted_spdmm_half": pred.get("spdmm_half", ""),
                   "predicted_full_rank":  pred.get("full_rank", "")}
            pd.DataFrame([row]).to_csv(csv_path, mode="a", header=False,
                                       index=False)
            print(f"  n={n:2d}  K={K:>3d}  method={info['method']:12s}  "
                  f"u={u_count:>6d}  cx={cx_count:>6d}", flush=True)


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="SIABLE benchmark: full-rank + low-rank dispatch + timing."
    )
    ap.add_argument("--nmin", type=int, default=3,
                    help="Minimum total qubit count (ancilla + data). Default: 3.")
    ap.add_argument("--nmax", type=int, default=7,
                    help="Maximum total qubit count. Default: 7. "
                         "Setting larger values is expensive (circuit size "
                         "grows as O(4^n)).")
    ap.add_argument("--reps", type=int, default=5,
                    help="Reps per full-rank compile-time cell.")
    ap.add_argument("--K-values", type=str, default="",
                    help="Comma-separated K values for the low-rank sweep. "
                         "Default: 1, 2, 4, 8, 16, 32, 64.")
    ap.add_argument("--skip-full-rank", action="store_true")
    ap.add_argument("--skip-low-rank",  action="store_true")
    ap.add_argument("--fresh", action="store_true",
                    help="Delete existing CSVs before running.")
    ap.add_argument("--append", action="store_true",
                    help="Append to existing CSVs (the default refuses to "
                         "overwrite pre-computed data).")
    args = ap.parse_args()

    n_values = list(range(args.nmin, args.nmax + 1))
    if args.K_values:
        K_values = [int(s) for s in args.K_values.split(",")]
    else:
        K_values = [1, 2, 4, 8, 16, 32, 64]

    if args.fresh:
        for p in (FULL_CSV, LOW_CSV):
            if os.path.exists(p):
                os.remove(p)
                print(f"  removed {p}")
    elif not args.append:
        existing = []
        if not args.skip_full_rank and os.path.exists(FULL_CSV):
            existing.append(FULL_CSV)
        if not args.skip_low_rank and os.path.exists(LOW_CSV):
            existing.append(LOW_CSV)
        if existing:
            print("[error] the following CSVs already exist:",
                  file=sys.stderr)
            for p in existing:
                print(f"          {p}", file=sys.stderr)
            print("        Re-run with --fresh to overwrite them, or with "
                  "--append to add rows.", file=sys.stderr)
            sys.exit(3)

    print(f"[config]  n = {n_values[0]}..{n_values[-1]}")
    print(f"[config]  K values for low-rank: {K_values}")
    print(f"[config]  reps (full-rank timing): {args.reps}")
    print(f"[config]  output CSVs: {FULL_CSV}, {LOW_CSV}")

    t0 = time.perf_counter()
    if not args.skip_full_rank:
        run_full_rank(n_values, args.reps, FULL_CSV)
    if not args.skip_low_rank:
        run_low_rank(n_values, K_values, LOW_CSV)
    elapsed = time.perf_counter() - t0

    print(f"\n[done] wall time: {elapsed:.1f} s ({elapsed / 60:.1f} min)")
    print(f"  full-rank -> {FULL_CSV}")
    print(f"  low-rank  -> {LOW_CSV}")


if __name__ == "__main__":
    main()
