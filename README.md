# SPDMM + SIABLE

Reference Python implementation of two algorithms from
**Li, Zhang & Zhang (2026)**,
*"Reducing C-NOT Counts for State Preparation and Block Encoding via Diagonal Matrix Migration"*
(arxiv `2603.16492`):

- **SPDMM** — *State Preparation via Diagonal Matrix Migration*.
  C-NOT leading constant `(11/12) · 2ⁿ`.
- **SIABLE** — *Single-Ancilla Block Encoding*.
  C-NOT leading constant `(11/48) · 4ⁿ⁺¹` (full rank).

The reference MATLAB / QCLAB code shipped with the paper is at
<https://github.com/zexianLIPolyU/SPDMM-SIABLE>; this repository ports
the same algorithms to Qiskit and provides the benchmark / figure
infrastructure used in the manuscript.

## Repository Layout

```
spdmm-siable/
├── spdmm/                            # the algorithm package
│   ├── __init__.py                   # public API
│   ├── spdmm.py                      # SPDMM state preparation
│   └── siable.py                     # SIABLE block encoding (full + low rank)
├── tests/                            # benchmarks + correctness tests
│   ├── test_siable.py                # pytest correctness suite (~31 s)
│   ├── _methods.py                   # state-prep method registry + helpers
│   ├── run_benchmark.py              # state-prep benchmark (gate counts + timing)
│   ├── run_siable_benchmark.py       # SIABLE benchmark (full + low rank + timing)
│   └── data/                         # pre-computed CSVs
│       ├── gate_counts.csv
│       ├── compile_time.csv
│       ├── siable_full_rank.csv
│       └── siable_low_rank.csv
├── plot_figure/                      # all figure-generation code
│   ├── _style.py                     # shared visual identity (palette, markers)
│   ├── plot_pareto.py                # state-prep Pareto frontier (n = 5, 10, 15)
│   ├── plot_compile_time.py          # state-prep compile-time scaling
│   └── plot_siable.py                # SIABLE scaling, dispatch, compile time
├── examples/                         # quick-start examples
│   ├── 01_state_prep.py
│   ├── 02_block_encoding.py
│   └── 03_low_rank.py
├── docs/                             # paper + supplementary LaTeX
├── pyproject.toml
├── requirements.txt
├── LICENSE                           # MIT
└── README.md
```

## Installation

```bash
pip install -e .[all]
```

The package itself only depends on Qiskit (>= 2.0), NumPy, SciPy, and
Numba.  The benchmark scripts additionally require `qclib` (for the
baselines), `pandas`, `matplotlib`, and `pytest`.

## Public API

```python
from spdmm import (
    state_preparation,        # n-qubit state preparation via SPDMM
    siable,                   # full-rank single-ancilla block encoding
    siable_low_rank,          # low-rank block encoding (auto-dispatch)
    SiableBlockEncoding,      # qiskit Gate wrapper for SIABLE
    SpdmmInitialize,          # qclib-style Initialize adapter for benchmarks
)
```

## Quick Start

### State preparation

```python
import numpy as np
from spdmm import state_preparation

n = 6
rng = np.random.default_rng(0)
psi = rng.standard_normal(2 ** n) + 1j * rng.standard_normal(2 ** n)
psi /= np.linalg.norm(psi)

qc = state_preparation(psi)
print(qc.count_ops())
```

### Block encoding

```python
from spdmm import siable
A = 0.3 * (rng.standard_normal((8, 8)) + 1j * rng.standard_normal((8, 8)))
qc, alpha = siable(A)
# Block-encoding relation: A == alpha * U[0::2, 0::2]
```

### Low-rank block encoding with auto-dispatch

```python
from spdmm import siable_low_rank
# A is a 16x16 rank-2 matrix
qc, alpha, info = siable_low_rank(A, K=2, return_info=True)
print(info["method"], info["nCNOT"])      # selected method + CNOT count
print(info["predicted_cost"])             # cost predicted for each method
```

## Reproducibility

The CSVs in `tests/data/` are the **authoritative measurements** used to
generate the figures in the paper.  They were produced with the package
versions pinned in `requirements.txt` (qiskit 2.4.1, qclib 0.1.15,
numpy 2.4.4, scipy 1.17.1, numba 0.65.1, pandas 3.0.2,
matplotlib 3.10.8).

For the published figures, please use the pre-computed CSVs directly —
running the benchmark from scratch with different package versions
will produce different absolute compile-time values (the relative
ordering and the gate counts are stable, but the absolute timings are
not).  See the Troubleshooting section below.

To pin the exact environment used for the original measurements:

```bash
pip install -r requirements.txt
```

To install with the looser forward-compatibility bounds (any qiskit
2.1+ / qclib 0.1.15+):

```bash
pip install -e .[all]
```

## Troubleshooting

### "Compile times look slower than the recorded CSV"

This is almost always an environment-drift issue, not a code change.
Three notes:

1. **Affected methods.** `low-rank`, `PB`, and `low-rank+BlockZXZ` route
   through qclib's pure-Python QSD internals (`qclib.unitary.unitary`
   / `decompose_unitary`).  These paths are *very* sensitive to qiskit
   minor-version updates, because qiskit's `HighLevelSynthesis` pass
   schedule can change which decomposition rules fire during
   `transpile()`.  In our own testing, a qiskit 2.0 → 2.4 jump made
   these three methods ~3× slower at n = 13..15.
2. **Unaffected methods.** `Isometry` (native qiskit), `PB+BlockZXZ`,
   `spdmm`, and `spdmm (Rust)` are essentially unchanged across
   qiskit versions — they don't touch the qclib QSD passes.  If only
   the qclib paths got slower, this is the cause.
3. **Recommended action.** For the paper figures, use the pre-computed
   CSVs.  For a sanity check on a fresh environment, run
   `python tests/_methods.py` to confirm every method produces the
   same gate counts as the CSV (the gate counts are stable; the times
   are not).

### "HighLevelSynthesis is unable to synthesize 'svd'"

This is an intermittent failure inside qclib's `SVDInitialize` /
`LowRankInitialize` that surfaces with newer qiskit versions
(observed with qiskit 2.4+).  qclib emits a custom gate labelled
`"svd"` and qiskit's newer `HighLevelSynthesis` pass occasionally
fails to find a decomposition rule for it.

The benchmark runner correctly drops the failed reps and reports the
remaining ones (you'll see `(4 reps)` instead of `(5 reps)` on the
affected lines).  This is an upstream qclib ↔ qiskit issue and does
not affect SPDMM or the Rust hybrid.

Workarounds: (a) pin to the versions in `requirements.txt`, (b)
re-run only the affected method with more reps, or (c) skip the
affected method via `--methods "spdmm,spdmm (Rust),Isometry,..."`.

### "Numba JIT takes a long time on the first SPDMM call"

The first invocation of `state_preparation` triggers a one-time Numba
compilation that can take 10-15 s.  Subsequent calls are fast
(milliseconds).  Both benchmark runners include an explicit warm-up
pass at the smallest n to amortise this; the self-test in
`tests/_methods.py` does not, so its first row will be slow.

## Reproducing the Paper Figures

The benchmark pipeline is two-stage: `tests/` produces CSV data and
`plot_figure/` reads it.

```bash
# Stage 1: regenerate the data (slow; ~30 min for full sweep)
cd tests/
python run_benchmark.py            # state-prep gate counts + compile time
python run_siable_benchmark.py     # full-rank + low-rank SIABLE

# Stage 2: regenerate the figures (fast; few seconds)
cd ../plot_figure/
python plot_pareto.py              # pareto_n5/n10/n15 + pareto_combined
python plot_compile_time.py        # compile_time
python plot_siable.py              # siable_fullrank_scaling + dispatch + compile_time
```

Pre-computed CSVs are checked into `tests/data/`, so step 2 works
without rerunning step 1.

## Methods Compared in the State-Prep Benchmark

| Method                | Reference                                    |
|-----------------------|----------------------------------------------|
| low-rank              | Low-rank state preparation (LRSP)            |
| PB                    | Plesch & Bužek                               |
| Isometry              | Iten et al. 2016                             |
| low-rank + BlockZXZ   | LRSP with Block-ZXZ for the inner unitaries  |
| PB + BlockZXZ         | PB with Block-ZXZ for the inner unitaries    |
| **spdmm**             | SPDMM (this work, pure-Python reference)     |
| spdmm (Rust)          | SPDMM hybrid using Qiskit's compiled kernel  |

The SPDMM family is rendered in bright deep pink across every figure
with the two variants distinguished by marker shape (star for the
pure-Python reference, filled plus for the Rust hybrid).

### The Rust hybrid

The `spdmm (Rust)` hybrid (described in §V-A of the paper) is bundled
as `spdmm/spdmm_hybrid_rust.py`.  It keeps the SPDMM outer recursion
and the Phase-3 isometry shell (Section III-B) intact, but replaces
the three leaf-level unitary syntheses inside `_isometry` and the
Phase-4 unitary with qiskit's compiled `qs_decomposition`.  The
diagonal-matrix-migration step is therefore lost: every Δ is identity.

Cost relative to pure-Python SPDMM at n = 15: **+7 CNOTs** (29634 vs
29627), **~2× faster** compile time.

The benchmark runner auto-detects the module at startup and includes
the hybrid as the seventh method.  If you want to swap in your own
hybrid implementation, replace this file with one exposing a class
`SpdmmInitializeRust` with the same `(state_vector,
label="spdmm (Rust)")` constructor and a `.definition` attribute
returning a QuantumCircuit.

## Safe Re-Run Workflow

By default, both benchmark runners **refuse to overwrite existing CSVs**
in `tests/data/`.  This protects the pre-computed paper data from being
silently polluted by partial runs.  Use one of:

- `--fresh` to delete the existing CSVs and start over from scratch
- `--append` to add rows to the existing CSVs (e.g. to extend the n
  range or add a new method)

## Running the Correctness Tests

```bash
pytest tests/test_siable.py -v        # 8 tests, ~30 s
```

## Citing

```bibtex
@misc{li2026reducingcnotcountsstate,
      title={Reducing C-NOT Counts for State Preparation and Block Encoding via Diagonal Matrix Migration}, 
      author={Zexian Li and Guofeng Zhang and Xiao-Ming Zhang},
      year={2026},
      eprint={2603.16492},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2603.16492}, 
}
```

## License

MIT (see `LICENSE`).
