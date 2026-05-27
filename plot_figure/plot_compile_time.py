"""
plot_compile_time.py
====================

Total compile time (build + transpile) versus qubit count n for the
seven state-preparation methods, using the visual identity defined in
``_style.py``.

Output (PNG at 300 dpi + vector PDF):
    compile_time.{png,pdf}

The data is read from ``../tests/data/compile_time.csv`` produced by
``tests/run_benchmark.py``.
"""

from __future__ import annotations

import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from _style import (
    METHODS_ORDER, STYLE_SINGLE,
    color_of, marker_of, is_spdmm,
)


_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "..", "tests", "data", "compile_time.csv")


def load_data(csv_path=_DEFAULT_CSV):
    """Return ``DATA[method][n] = (mean_ms, ci_ms)``."""
    if not os.path.exists(csv_path):
        raise SystemExit(
            f"Compile-time CSV not found at:\n  {csv_path}\n"
            "Run `python ../tests/run_benchmark.py` first.")
    df = pd.read_csv(csv_path)
    out = {}
    for m in df["method"].unique():
        sub = df[df["method"] == m].set_index("n").sort_index()
        out[m] = {
            int(n): (float(sub.loc[n, "total_mean"]) * 1000.0,
                     float(sub.loc[n, "total_ci95"]) * 1000.0)
            for n in sub.index
        }
    return out


# =============================================================================
# Figure
# =============================================================================

def make_figure(out_path="compile_time"):
    DATA = load_data()
    plt.rcParams.update(STYLE_SINGLE)
    plt.rcParams.update({"axes.labelsize": 18, "xtick.labelsize": 15,
                         "ytick.labelsize": 15, "legend.fontsize": 13})

    fig, ax = plt.subplots(figsize=(7.8, 5.0))

    # SPDMM halo discs drawn first so markers sit on top
    for m in METHODS_ORDER:
        if m not in DATA or not is_spdmm(m):
            continue
        rec = sorted(DATA[m].items())
        ns = [r[0] for r in rec]
        means = [r[1][0] for r in rec]
        ax.scatter(ns, means, marker="o", s=420,
                   facecolor=color_of(m), edgecolor="none",
                   alpha=0.16, zorder=2)

    # Line + error bars for every method
    for m in METHODS_ORDER:
        if m not in DATA:
            continue
        rec = sorted(DATA[m].items())
        ns = [r[0] for r in rec]
        means = [r[1][0] for r in rec]
        cis = [r[1][1] for r in rec]
        is_sp = is_spdmm(m)
        ax.errorbar(
            ns, means, yerr=cis,
            fmt=marker_of(m) + "-",
            color=color_of(m),
            capsize=4, capthick=1.8 if is_sp else 1.5,
            lw=3.0 if is_sp else 2.0,
            markersize=16 if is_sp else 9,
            markeredgewidth=1.8 if is_sp else 0.6,
            markeredgecolor="black" if is_sp else "white",
            zorder=4 if is_sp else 3,
            label=m,
        )

    ax.set_xlabel(r"number of qubits $n$")
    ax.set_ylabel("total compile time (ms, log scale)")
    ax.set_yscale("log")
    ax.grid(True, which="both", ls=":", alpha=0.55)
    ax.set_xticks(sorted(set(n for d in DATA.values() for n in d)))

    # Custom legend mirroring the in-plot emphasis convention
    handles = []
    for m in METHODS_ORDER:
        if m not in DATA:
            continue
        is_sp = is_spdmm(m)
        handles.append(Line2D(
            [0], [0],
            marker=marker_of(m),
            color=color_of(m),
            markerfacecolor=color_of(m),
            markeredgecolor="black" if is_sp else "white",
            markeredgewidth=1.5 if is_sp else 0.5,
            markersize=15 if is_sp else 10,
            lw=2.8 if is_sp else 2.0,
            label=m,
        ))
    ax.legend(handles=handles, framealpha=0.95, loc="upper left", ncol=2,
              handlelength=2.4, columnspacing=1.2, borderaxespad=0.4)

    # "Better" arrow pointing toward the SPDMM curves
    ax.annotate("", xy=(0.95, 0.05), xytext=(0.82, 0.18),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=2.2))
    ax.text(0.95, 0.20, "better", transform=ax.transAxes,
            fontsize=11, color="#2ca02c", fontweight="bold", ha="right")

    fig.tight_layout()
    fig.savefig(f"{out_path}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_path}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_path}.png / .pdf")


if __name__ == "__main__":
    make_figure()
