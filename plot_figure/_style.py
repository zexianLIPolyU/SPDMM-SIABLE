"""
Shared visual identity for all benchmark figures.

All plotting scripts in this folder import constants from this module
so the same method receives the same colour, marker, and font in every
figure across the paper.
"""

# -----------------------------------------------------------------------------
# State-preparation methods
# -----------------------------------------------------------------------------

PALETTE = {
    # Baselines (Seaborn tab10 muted palette)
    "low-rank":          "#4C72B0",   # blue
    "PB":                "#DD8452",   # orange
    "Isometry":          "#55A868",   # green
    "low-rank+BlockZXZ": "#8172B3",   # purple
    "PB+BlockZXZ":       "#C44E52",   # red
    # SPDMM family: shared bright pink, distinguished by marker
    "spdmm":             "#FF1493",
    "spdmm (Rust)":      "#FF1493",
}

MARKER = {
    "low-rank":          "o",
    "PB":                "s",
    "Isometry":          "^",
    "low-rank+BlockZXZ": "D",
    "PB+BlockZXZ":       "v",
    "spdmm":             "*",
    "spdmm (Rust)":      "P",
}

METHODS_ORDER = [
    "low-rank", "PB", "Isometry",
    "low-rank+BlockZXZ", "PB+BlockZXZ",
    "spdmm", "spdmm (Rust)",
]

# Subset used by the Pareto figures (no Rust hybrid)
METHODS_PARETO = [m for m in METHODS_ORDER if m != "spdmm (Rust)"]


# -----------------------------------------------------------------------------
# SIABLE methods (low-rank dispatch)
# -----------------------------------------------------------------------------
# The four sub-constructions selected by ``siable_low_rank(method='auto')``.
# Colour scheme intentionally echoes the state-prep palette so a reader
# recognises ``state_prep`` (uses SPDMM internally) as pink and
# ``isometry`` (uses qiskit Isometry) as green.

SIABLE_PALETTE = {
    "state_prep": "#FF1493",   # bright pink (uses SPDMM)
    "isometry":   "#55A868",   # green (uses qiskit Isometry)
    "spdmm_half": "#8172B3",   # purple (Block-ZXZ half-isometry hybrid)
    "full_rank":  "#DD8452",   # orange (general Block-ZXZ unitary path)
}
SIABLE_MARKER = {
    "state_prep": "*",
    "isometry":   "^",
    "spdmm_half": "D",
    "full_rank":  "s",
}

# Pretty display names for figures and tables
SIABLE_LABEL = {
    "state_prep": "state-prep ($K=1$)",
    "isometry":   "isometry",
    "spdmm_half": "SPDMM half-iso.",
    "full_rank":  "full-rank",
}
SIABLE_METHODS_ORDER = ["state_prep", "isometry", "spdmm_half", "full_rank"]


# -----------------------------------------------------------------------------
# Matplotlib rcParams
# -----------------------------------------------------------------------------

STYLE_SINGLE = {
    "font.family":       "serif",
    "font.size":         14,
    "axes.labelsize":    17,
    "xtick.labelsize":   13,
    "ytick.labelsize":   13,
    "legend.fontsize":   11,
    "axes.linewidth":    1.3,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size":  5,
    "ytick.major.size":  5,
}

STYLE_MULTI = {
    "font.family":       "serif",
    "font.size":         13,
    "axes.labelsize":    15,
    "xtick.labelsize":   12,
    "ytick.labelsize":   12,
    "legend.fontsize":   12,
    "axes.linewidth":    1.2,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "xtick.major.size":  4,
    "ytick.major.size":  4,
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def color_of(method):
    """Look up the canonical colour for any state-prep or SIABLE method."""
    return PALETTE.get(method) or SIABLE_PALETTE.get(method) or "#000000"


def marker_of(method):
    """Look up the canonical marker shape for any method."""
    return MARKER.get(method) or SIABLE_MARKER.get(method) or "o"


def is_spdmm(method):
    """True if the method belongs to the SPDMM family (state_prep included)."""
    m = method.lower()
    return m.startswith("spdmm") or m == "state_prep"


# ---------------------------------------------------------------------------
# Self-test (run this file directly to inspect the palette)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("State-preparation methods")
    print(f"  {'method':<22s} {'colour':<10s} {'marker':<8s}")
    for m in METHODS_ORDER:
        print(f"  {m:<22s} {PALETTE[m]:<10s} {MARKER[m]:<8s}")
    print()
    print("SIABLE low-rank sub-methods (dispatch)")
    print(f"  {'method':<14s} {'colour':<10s} {'marker':<8s}")
    for m in SIABLE_METHODS_ORDER:
        print(f"  {m:<14s} {SIABLE_PALETTE[m]:<10s} {SIABLE_MARKER[m]:<8s}")
