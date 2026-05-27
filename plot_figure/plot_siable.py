"""
plot_siable.py
==============

Figures for the SIABLE (Single Ancilla Block Encoding) experiments.
Reads from CSVs produced by ``tests/run_siable_benchmark.py``.

Three figures are generated:

  1. ``siable_fullrank_scaling.{png,pdf}``
        Full-rank SIABLE CNOT count vs n, with the theoretical leading
        constant (11/12) * 2^(n+1) annotated.

  2. ``siable_low_rank_dispatch.{png,pdf}``
        At the largest available n, CNOT count vs rank K for each of the
        four sub-constructions (state_prep / isometry / spdmm_half /
        full_rank), with the auto-dispatch envelope highlighted.

  3. ``siable_compile_time.{png,pdf}``
        Full-rank SIABLE compile time vs n with 95% CI error bars.

All three figures share the visual identity defined in ``_style.py``
so the SIABLE colour-coded methods match the same methods that appear
on the state-preparation Pareto / compile-time figures (state_prep is
pink, isometry is green, etc.).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from _style import (
    SIABLE_PALETTE, SIABLE_MARKER, SIABLE_LABEL, SIABLE_METHODS_ORDER,
    STYLE_SINGLE,
)


_HERE = os.path.dirname(os.path.abspath(__file__))
_FULL_CSV = os.path.join(_HERE, "..", "tests", "data", "siable_full_rank.csv")
_LOW_CSV  = os.path.join(_HERE, "..", "tests", "data", "siable_low_rank.csv")


def _load_full_rank():
    if not os.path.exists(_FULL_CSV):
        raise SystemExit(
            f"Full-rank SIABLE CSV not found at:\n  {_FULL_CSV}\n"
            "Run `python ../tests/run_siable_benchmark.py` first.")
    return pd.read_csv(_FULL_CSV).sort_values("n")


def _load_low_rank():
    if not os.path.exists(_LOW_CSV):
        raise SystemExit(
            f"Low-rank SIABLE CSV not found at:\n  {_LOW_CSV}\n"
            "Run `python ../tests/run_siable_benchmark.py` first.")
    return pd.read_csv(_LOW_CSV).sort_values(["n", "K"])


# =============================================================================
# Figure 1: Full-rank SIABLE scaling
# =============================================================================

def fig_full_rank_scaling(out_path="siable_fullrank_scaling"):
    df = _load_full_rank()
    plt.rcParams.update(STYLE_SINGLE)
    plt.rcParams.update({"axes.labelsize": 17, "legend.fontsize": 12})

    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    ns = df["n"].to_numpy()
    cx = df["cx"].to_numpy()
    u  = df["u"].to_numpy()

    # CNOT curve
    ax.plot(ns, cx, "P-", color="#FF1493",
            lw=2.8, markersize=14, markeredgewidth=1.4,
            markeredgecolor="black", zorder=4,
            label=r"SIABLE CNOT count (measured)")

    # Single-qubit curve, plotted as a thinner ghost for context
    ax.plot(ns, u, "o-", color="#888888", lw=1.6, markersize=8,
            alpha=0.85, label=r"SIABLE single-qubit count")

    # Theoretical leading constant: (11/12) * 2^(n+1) = (11/48) * 4^(n+1) ?
    # Paper Table III leading constant is (11/48) * 4^(n+1) total CNOTs.
    # That equals (11/12) * 4^n in terms of n.  Compute on a smooth grid:
    ns_smooth = np.linspace(ns.min(), ns.max(), 100)
    theory = (11.0 / 48.0) * 4 ** (ns_smooth + 1)
    ax.plot(ns_smooth, theory, "--", color="#444444", lw=1.5,
            label=r"theoretical leading constant $\frac{11}{48}\cdot4^{n+1}$")

    ax.set_xlabel(r"total qubit count $n$ (1 ancilla + $n{-}1$ data)")
    ax.set_ylabel("gate count")
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.grid(True, which="both", ls=":", alpha=0.55)
    ax.legend(loc="upper left", framealpha=0.95)

    # Annotate the measured CNOT counts inline
    for ni, ci in zip(ns, cx):
        ax.annotate(f"{ci}", (ni, ci), xytext=(0, 12),
                    textcoords="offset points",
                    fontsize=10, ha="center", color="#C71585")

    fig.tight_layout()
    fig.savefig(f"{out_path}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_path}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_path}.png / .pdf")


# =============================================================================
# Figure 2: Low-rank SIABLE method dispatch at fixed n
# =============================================================================

def fig_low_rank_dispatch(out_path="siable_low_rank_dispatch", n=None):
    df = _load_low_rank()
    if n is None:
        n = int(df["n"].max())
    sub = df[df["n"] == n].sort_values("K").reset_index(drop=True)
    if sub.empty:
        print(f"  [warn] no low-rank data at n={n}; skipping figure")
        return

    plt.rcParams.update(STYLE_SINGLE)
    plt.rcParams.update({"axes.labelsize": 17, "legend.fontsize": 12})
    fig, ax = plt.subplots(figsize=(8.0, 5.0))

    Ks = sub["K"].to_numpy()

    # For each method, plot its predicted cost across all K values
    # (these come from info['predicted_cost'] inside the benchmark runner)
    for m in SIABLE_METHODS_ORDER:
        col = f"predicted_{m}"
        if col not in sub.columns:
            continue
        ys = pd.to_numeric(sub[col], errors="coerce").to_numpy()
        finite = np.isfinite(ys)
        if finite.sum() == 0:
            continue
        ax.plot(Ks[finite], ys[finite],
                marker=SIABLE_MARKER[m], color=SIABLE_PALETTE[m],
                lw=1.8, markersize=10, alpha=0.7,
                label=f"predicted: {SIABLE_LABEL[m]}")

    # Overlay the actually-selected method per K with bold, larger markers
    for _, row in sub.iterrows():
        K = row["K"]; cx = int(row["cx"]); chosen = row["method"]
        ax.scatter([K], [cx], marker=SIABLE_MARKER[chosen],
                   color=SIABLE_PALETTE[chosen],
                   s=400 if chosen == "state_prep" else 260,
                   edgecolors="black", linewidths=2.0, zorder=5)

    # Build a legend that also explains the "selected by auto" marker
    chosen_handles = []
    for m in SIABLE_METHODS_ORDER:
        chosen_handles.append(Line2D(
            [0], [0], marker=SIABLE_MARKER[m], color="white",
            markerfacecolor=SIABLE_PALETTE[m],
            markeredgecolor="black", markeredgewidth=1.6,
            markersize=14, label=f"selected: {SIABLE_LABEL[m]}",
        ))
    leg1 = ax.legend(loc="upper left", framealpha=0.95, fontsize=10,
                     title=r"predicted CNOT (per method)", title_fontsize=10)
    ax.add_artist(leg1)
    ax.legend(handles=chosen_handles, loc="lower right",
              framealpha=0.95, fontsize=10,
              title=r"auto-dispatch choice", title_fontsize=10)

    ax.set_xlabel(r"rank $K$")
    ax.set_ylabel("CNOT count")
    ax.set_yscale("log")
    ax.set_xscale("log", base=2)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.set_title(f"SIABLE low-rank method dispatch at $n = {n}$",
                 fontsize=14, pad=8)

    fig.tight_layout()
    fig.savefig(f"{out_path}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_path}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_path}.png / .pdf")


# =============================================================================
# Figure 3: SIABLE compile time
# =============================================================================

def fig_compile_time(out_path="siable_compile_time"):
    df = _load_full_rank()
    plt.rcParams.update(STYLE_SINGLE)
    plt.rcParams.update({"axes.labelsize": 17, "legend.fontsize": 12})
    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    ns = df["n"].to_numpy()
    total_ms = df["total_mean"].to_numpy() * 1000.0
    ci_ms = df["total_ci95"].to_numpy() * 1000.0

    # Halo
    ax.scatter(ns, total_ms, marker="o", s=420,
               facecolor="#FF1493", edgecolor="none", alpha=0.16, zorder=2)
    # Main curve
    ax.errorbar(ns, total_ms, yerr=ci_ms,
                fmt="P-", color="#FF1493",
                capsize=4, capthick=1.8,
                lw=3.0, markersize=16,
                markeredgewidth=1.8, markeredgecolor="black",
                zorder=4, label="SIABLE compile time")

    ax.set_xlabel(r"total qubit count $n$")
    ax.set_ylabel("compile time (ms, log scale)")
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.grid(True, which="both", ls=":", alpha=0.55)
    ax.legend(loc="upper left", framealpha=0.95)

    fig.tight_layout()
    fig.savefig(f"{out_path}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_path}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_path}.png / .pdf")


if __name__ == "__main__":
    fig_full_rank_scaling()
    fig_low_rank_dispatch()
    fig_compile_time()
