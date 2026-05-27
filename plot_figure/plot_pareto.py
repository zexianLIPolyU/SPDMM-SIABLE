"""
plot_pareto.py
==============

Pareto frontier figures in the (single-qubit count, CNOT count) plane
at qubit counts n = 5, 10, and 15.  Only the six primary state-prep
methods are included; the SPDMM-Rust hybrid is omitted because it is a
software engineering convenience rather than a distinct algorithm and
would visually duplicate the ``spdmm`` point on these plots.

Outputs:
    pareto_n5.{png,pdf}
    pareto_n10.{png,pdf}
    pareto_n15.{png,pdf}
    pareto_combined.{png,pdf}     -- three-panel composite with shared legend

The data is read from ``../tests/data/gate_counts.csv`` produced by
``tests/run_benchmark.py``.
"""

from __future__ import annotations

import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

from _style import (
    METHODS_PARETO, STYLE_SINGLE, STYLE_MULTI,
    color_of, marker_of, is_spdmm,
)


_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CSV = os.path.join(_HERE, "..", "tests", "data", "gate_counts.csv")


def load_data(csv_path=_DEFAULT_CSV):
    """Return ``DATA[method] = {'u': {n: count}, 'cx': {n: count}}``."""
    if not os.path.exists(csv_path):
        raise SystemExit(
            f"Gate-count CSV not found at:\n  {csv_path}\n"
            "Run `python ../tests/run_benchmark.py` first.")
    df = pd.read_csv(csv_path)
    data = {}
    for m in df["method"].unique():
        sub = df[df["method"] == m].set_index("n").sort_index()
        data[m] = {
            "u":  {int(n): int(sub.loc[n, "u"])  for n in sub.index},
            "cx": {int(n): int(sub.loc[n, "cx"]) for n in sub.index},
        }
    return data


LABEL_OFFSETS = {
    5: {
        "spdmm":             (18,   22),
        "Isometry":          (16,   10),
        "low-rank":          (16,   10),
        "low-rank+BlockZXZ": (-92, -36),
        "PB+BlockZXZ":       (-105, 22),
        "PB":                (16,   10),
    },
    10: {
        "spdmm":             (-22, -50),
        "low-rank+BlockZXZ": (65,   30),
        "PB+BlockZXZ":       (-115, 28),
        "low-rank":          (-90,  35),
        "PB":                (40,  -38),
        "Isometry":          (16,   10),
    },
    15: {
        "spdmm":             (16,   22),
        "Isometry":          (16,   10),
        "low-rank":          (-92, -22),
        "low-rank+BlockZXZ": (16,   28),
        "PB+BlockZXZ":       (-118, 18),
        "PB":                (-30, -28),
    },
}


def pareto_indices(points):
    n = len(points)
    dominated = [False] * n
    for i in range(n):
        xi, yi = points[i]
        for j in range(n):
            if i == j:
                continue
            xj, yj = points[j]
            if xj <= xi and yj <= yi and (xj < xi or yj < yi):
                dominated[i] = True
                break
    optimal = [i for i in range(n) if not dominated[i]]
    optimal.sort(key=lambda i: points[i][0])
    return optimal


def plot_panel(ax, DATA, n, show_labels=True):
    x_vals = [DATA[m]["u"][n]  for m in METHODS_PARETO]
    y_vals = [DATA[m]["cx"][n] for m in METHODS_PARETO]
    points = list(zip(x_vals, y_vals))
    pareto = set(pareto_indices(points))

    x_pad = (max(x_vals) - min(x_vals)) * 0.16
    y_pad = (max(y_vals) - min(y_vals)) * 0.16
    xlim = (min(x_vals) - x_pad, max(x_vals) + x_pad * 1.6)
    ylim = (min(y_vals) - y_pad, max(y_vals) + y_pad)
    ax.set_xlim(xlim); ax.set_ylim(ylim)

    # Dominated region (union of upper-right quadrants from each Pareto point)
    pf_points = sorted([points[i] for i in pareto])
    if len(pf_points) >= 2:
        poly = [(pf_points[0][0], ylim[1])]
        poly.append((xlim[1], ylim[1]))
        poly.append((xlim[1], pf_points[-1][1]))
        for i in range(len(pf_points) - 1, -1, -1):
            poly.append(pf_points[i])
            if i > 0:
                poly.append((pf_points[i][0], pf_points[i - 1][1]))
        ax.add_patch(Polygon(poly, closed=True,
                             facecolor="#d62728", alpha=0.05,
                             edgecolor="none", zorder=1))
        for i in range(len(pf_points) - 1):
            x1, y1 = pf_points[i]; x2, y2 = pf_points[i + 1]
            ax.plot([x1, x2], [y1, y1], "--", color="#444444",
                    lw=2.0, zorder=2)
            ax.plot([x2, x2], [y1, y2], "--", color="#444444",
                    lw=2.0, zorder=2)
        ax.plot([xlim[0], pf_points[0][0]], [pf_points[0][1]] * 2,
                ":", color="#888888", lw=1.2, zorder=2)
        ax.plot([pf_points[-1][0]] * 2, [pf_points[-1][1], ylim[1]],
                ":", color="#888888", lw=1.2, zorder=2)

    # SPDMM emphasis: double halo
    spdmm_idx = METHODS_PARETO.index("spdmm")
    sx, sy = x_vals[spdmm_idx], y_vals[spdmm_idx]
    ax.scatter(sx, sy, marker="o", s=1100,
               facecolor="#FFD700", edgecolor="none", alpha=0.20, zorder=2.5)
    ax.scatter(sx, sy, marker="o", s=650,
               facecolor=color_of("spdmm"), edgecolor="none",
               alpha=0.30, zorder=2.8)

    # Data points
    for i, m in enumerate(METHODS_PARETO):
        x, y = x_vals[i], y_vals[i]
        is_opt = (i in pareto)
        is_sp = is_spdmm(m)
        if is_sp:
            ax.scatter(x, y, marker=marker_of(m), color=color_of(m),
                       s=500, edgecolors="black", linewidths=2.5, zorder=5)
        elif is_opt:
            ax.scatter(x, y, marker=marker_of(m), color=color_of(m),
                       s=260, edgecolors="black", linewidths=1.8, zorder=5)
        else:
            ax.scatter(x, y, marker=marker_of(m), color=color_of(m),
                       s=140, alpha=0.6,
                       edgecolors="white", linewidths=0.6, zorder=4)

    # Labels with leader lines
    if show_labels:
        offsets = LABEL_OFFSETS[n]
        for i, m in enumerate(METHODS_PARETO):
            x, y = x_vals[i], y_vals[i]
            dx, dy = offsets.get(m, (12, 8))
            is_sp = is_spdmm(m)
            ax.annotate(
                m, xy=(x, y), xytext=(dx, dy),
                textcoords="offset points",
                fontsize=13 if is_sp else 11,
                fontweight="bold" if (is_sp or m == "Isometry") else "normal",
                color=color_of(m),
                bbox=dict(
                    boxstyle="round,pad=0.28" if is_sp else "round,pad=0.18",
                    facecolor="#FFF59D" if is_sp else "white",
                    edgecolor="#F57F17" if is_sp else "#bbbbbb",
                    linewidth=1.6 if is_sp else 0.6, alpha=0.95,
                ),
                arrowprops=dict(arrowstyle="-", color="#888888",
                                lw=0.7, shrinkA=4, shrinkB=6),
                zorder=6,
            )

    ax.annotate("", xy=(0.05, 0.05), xytext=(0.18, 0.18),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=2.2))
    ax.text(0.05, 0.20, "better", transform=ax.transAxes,
            fontsize=11, color="#2ca02c", fontweight="bold")

    ax.set_xlabel(r"single-qubit ($u$) gate count")
    ax.set_ylabel("CNOT (cx) gate count")
    ax.grid(True, ls=":", alpha=0.45, zorder=0)

    return pareto


def add_panel_legend(ax, pareto):
    handles = []
    for i, m in enumerate(METHODS_PARETO):
        is_opt = (i in pareto)
        is_sp = is_spdmm(m)
        handles.append(Line2D(
            [0], [0], marker=marker_of(m), color="white",
            markerfacecolor=color_of(m),
            markeredgecolor="black" if (is_opt or is_sp) else "white",
            markeredgewidth=2.0 if is_sp else (1.5 if is_opt else 0.5),
            markersize=17 if is_sp else (15 if is_opt else 10),
            label=m,
        ))
    handles.append(Line2D([0], [0], linestyle="--", color="#444444",
                          lw=2.0, label="Pareto frontier"))
    ax.legend(handles=handles, loc="upper right",
              framealpha=0.95, fontsize=11, ncol=1, borderaxespad=0.4)


def make_single_panel(n, DATA, out_prefix="pareto"):
    plt.rcParams.update(STYLE_SINGLE)
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    pareto = plot_panel(ax, DATA, n)
    add_panel_legend(ax, pareto)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_n{n}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_n{n}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_prefix}_n{n}.png / .pdf")


def make_combined(DATA, ns=(5, 10, 15), out_prefix="pareto"):
    plt.rcParams.update(STYLE_MULTI)
    fig, axes = plt.subplots(1, len(ns), figsize=(6.0 * len(ns), 5.6))
    for ax, n in zip(axes, ns):
        plot_panel(ax, DATA, n, show_labels=False)
        ax.set_title(f"$n = {n}$", fontsize=15, pad=8)

    optimal_any = set()
    for n in ns:
        xs = [DATA[m]["u"][n]  for m in METHODS_PARETO]
        ys = [DATA[m]["cx"][n] for m in METHODS_PARETO]
        for idx in pareto_indices(list(zip(xs, ys))):
            optimal_any.add(METHODS_PARETO[idx])

    handles = []
    for m in METHODS_PARETO:
        is_opt = m in optimal_any
        is_sp = is_spdmm(m)
        handles.append(Line2D(
            [0], [0], marker=marker_of(m), color="white",
            markerfacecolor=color_of(m),
            markeredgecolor="black" if (is_opt or is_sp) else "white",
            markeredgewidth=2.2 if is_sp else (1.5 if is_opt else 0.5),
            markersize=18 if is_sp else (15 if is_opt else 10),
            label=m,
        ))
    handles.append(Line2D([0], [0], linestyle="--", color="#444444",
                          lw=2.0, label="Pareto frontier"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles),
               bbox_to_anchor=(0.5, 1.04), frameon=False, fontsize=12.5,
               handlelength=2.2, columnspacing=1.6)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{out_prefix}_combined.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_prefix}_combined.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_prefix}_combined.png / .pdf")


if __name__ == "__main__":
    DATA = load_data()
    for n in (5, 10, 15):
        make_single_panel(n, DATA)
    make_combined(DATA)
