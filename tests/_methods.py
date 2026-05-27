"""
Method registry and benchmark helpers for state-preparation experiments.

Imports the six "primary" state-preparation methods used in the paper
benchmarks; an optional ``spdmm (Rust)`` hybrid is picked up if the user
provides it at the import path ``spdmm.spdmm_hybrid_rust`` with a class
``SpdmmInitializeRust``.

Only consumed by ``run_benchmark.py``; the plotting scripts do not
import this module.
"""

from __future__ import annotations

import importlib
import os
import pkgutil
import sys
import time

import numpy as np
from scipy import stats as sstats

import qiskit
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import Isometry as _QiskitIsometry
from qiskit.synthesis.unitary.qsd import qs_decomposition

import qclib.unitary as _qu  # noqa: F401  (side-effects)
import qclib.state_preparation as _sp_pkg
from qclib.state_preparation import LowRankInitialize, SVDInitialize

# Add the repository root to sys.path so we can import the local spdmm package
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from spdmm import SpdmmInitialize  # noqa: E402


# ---------------------------------------------------------------------------
# Native-qiskit Isometry adapter (qclib-style)
# ---------------------------------------------------------------------------

class IsometryInitialize:
    """qclib-style adapter around qiskit's native ``Isometry`` gate.

    Matches the ``isometry`` baseline used in qclib-papers
    state-preparation benchmark notebook; all synthesis work is performed
    inside qiskit during transpile, which is ~100x faster than qclib's
    pure-Python ``scheme='ccd'`` mode.
    """

    def __init__(self, state_vector, label="iso"):
        v = np.asarray(state_vector, dtype=complex).reshape(-1)
        n = int(round(np.log2(len(v))))
        v = v / np.linalg.norm(v)
        iso = _QiskitIsometry(v.reshape(-1, 1), 0, 0)
        qc = QuantumCircuit(n)
        qc.append(iso, qc.qubits)
        self.definition = qc
        self.num_qubits = n


# ---------------------------------------------------------------------------
# Block-ZXZ monkey-patch infrastructure
# ---------------------------------------------------------------------------
# qclib's LowRankInitialize and SVDInitialize import ``unitary`` /
# ``decompose_unitary`` from ``qclib.unitary``.  Replacing those names
# with a thin wrapper around qiskit's ``qs_decomposition`` (= Block-ZXZ
# in qiskit >= 2.1) turns the methods into their "+BlockZXZ" variants.

if "Block ZXZ" not in (qs_decomposition.__doc__ or ""):
    raise RuntimeError(
        f"qiskit {qiskit.__version__} qs_decomposition is NOT Block-ZXZ. "
        "Please upgrade to qiskit >= 2.1.")
print(f"[OK] qiskit {qiskit.__version__}, qs_decomposition = Block-ZXZ")

_PATCH_NAMES = ("decompose_unitary", "unitary")


def _find_patch_targets():
    out = {}
    for info in pkgutil.iter_modules(_sp_pkg.__path__,
                                     prefix=_sp_pkg.__name__ + "."):
        try:
            mod = importlib.import_module(info.name)
        except Exception:
            continue
        for attr in _PATCH_NAMES:
            if hasattr(mod, attr):
                fn = getattr(mod, attr)
                if getattr(fn, "__module__", "") == "qclib.unitary":
                    out[(mod, attr)] = fn
    return out


_TARGETS = _find_patch_targets()
print(f"[OK] Found {len(_TARGETS)} qclib patch points "
      f"(low-rank / SVD will be re-routed in '+BlockZXZ' variants)")


def _block_zxz_unitary(gate, decomposition="qsd",
                       opt_a1=True, opt_a2=True, **_):
    """Replacement for ``qclib.unitary.unitary`` that uses qiskit's
    Block-ZXZ ``qs_decomposition`` instead of qclib's pure-Python KAK."""
    mat = np.asarray(gate, dtype=complex)
    n = int(np.log2(mat.shape[0]))
    qc = QuantumCircuit(n)
    sub = qs_decomposition(mat, opt_a1=opt_a1, opt_a2=opt_a2)
    qc.append(sub.to_gate(label="block_zxz"), qc.qubits)
    return qc.to_gate()


class UnitaryBackend:
    """Context manager that patches qclib to use Block-ZXZ when needed."""

    def __init__(self, backend):
        self.backend = backend

    def __enter__(self):
        if self.backend == "block_zxz":
            for (mod, attr), _orig in _TARGETS.items():
                setattr(
                    mod, attr,
                    lambda g, decomposition="qsd", **k:
                        _block_zxz_unitary(g, decomposition, **k),
                )
        else:
            for (mod, attr), orig in _TARGETS.items():
                setattr(mod, attr, orig)
        return self

    def __exit__(self, *exc):
        for (mod, attr), orig in _TARGETS.items():
            setattr(mod, attr, orig)


# ---------------------------------------------------------------------------
# Optional Rust hybrid
# ---------------------------------------------------------------------------

_SpdmmRust = None
try:
    _hybrid_mod = importlib.import_module("spdmm.spdmm_hybrid_rust")
    _SpdmmRust = getattr(_hybrid_mod, "SpdmmInitializeRust", None)
    if _SpdmmRust is not None:
        print("[OK] spdmm.spdmm_hybrid_rust found -- 'spdmm (Rust)' enabled")
except ImportError:
    pass

if _SpdmmRust is None:
    print("[note] spdmm.spdmm_hybrid_rust not found -- "
          "'spdmm (Rust)' will be skipped by the benchmark runner")


# ---------------------------------------------------------------------------
# Method registry
# ---------------------------------------------------------------------------

METHODS = {
    "low-rank":          (LowRankInitialize,  "qsd"),
    "PB":                (SVDInitialize,      "qsd"),
    "Isometry":          (IsometryInitialize, "qsd"),
    "low-rank+BlockZXZ": (LowRankInitialize,  "block_zxz"),
    "PB+BlockZXZ":       (SVDInitialize,      "block_zxz"),
    "spdmm":             (SpdmmInitialize,    "block_zxz"),
}
if _SpdmmRust is not None:
    METHODS["spdmm (Rust)"] = (_SpdmmRust, "block_zxz")


# ---------------------------------------------------------------------------
# Benchmark helpers
# ---------------------------------------------------------------------------

def random_state(n, seed):
    rng = np.random.RandomState(seed)
    v = rng.randn(2 ** n) + 1j * rng.randn(2 ** n)
    return v / np.linalg.norm(v)


def count_gates(n, seed, method):
    """Build the circuit and return (u_count, cx_count) after transpilation."""
    psi = random_state(n, seed)
    InitCls, backend = METHODS[method]
    with UnitaryBackend(backend):
        init = InitCls(psi)
        circ = init.definition
    tc = transpile(circ, basis_gates=["u", "cx"], optimization_level=0)
    ops = tc.count_ops()
    return int(ops.get("u", 0)), int(ops.get("cx", 0))


def time_one(n, seed, method):
    """Return (t_build, t_transpile, cx_count) for one cell."""
    psi = random_state(n, seed)
    InitCls, backend = METHODS[method]
    with UnitaryBackend(backend):
        t0 = time.perf_counter()
        init = InitCls(psi)
        circ = init.definition
        t_build = time.perf_counter() - t0
    t0 = time.perf_counter()
    tc = transpile(circ, basis_gates=["u", "cx"], optimization_level=0)
    t_transpile = time.perf_counter() - t0
    return t_build, t_transpile, tc.count_ops().get("cx", 0)


def mean_ci(samples, alpha=0.05):
    a = np.asarray(samples, dtype=float)
    n = len(a)
    m = a.mean()
    if n < 2:
        return m, 0.0
    se = a.std(ddof=1) / np.sqrt(n)
    h = se * sstats.t.ppf(1 - alpha / 2, n - 1)
    return m, h


# ---------------------------------------------------------------------------
# Self-test (run this file directly to verify all methods are functional)
# ---------------------------------------------------------------------------

def _self_test(n=4, seed=0):
    """Build the circuit for each method in METHODS once at the given n,
    transpile to the {u, cx} basis, and report counts + timing.

    Useful as a quick sanity check after pulling the repo or after
    adding/replacing the Rust hybrid -- it touches every method in the
    registry and surfaces any setup errors immediately.
    """
    print()
    print(f"=== self-test: n={n}, seed={seed} ===")
    print(f"  {'method':<20s} {'u':>7s} {'cx':>7s} {'build(ms)':>11s} "
          f"{'transpile(ms)':>14s}")
    failures = []
    for method in METHODS:
        try:
            tb, tt, cx = time_one(n, seed, method)
            u, cx2 = count_gates(n, seed, method)
            assert cx == cx2, f"timing/count mismatch: {cx} vs {cx2}"
            print(f"  {method:<20s} {u:>7d} {cx:>7d} "
                  f"{tb * 1000:>11.2f} {tt * 1000:>14.2f}")
        except Exception as e:
            print(f"  {method:<20s} FAIL: {type(e).__name__}: {e}")
            failures.append(method)
    print()
    if failures:
        print(f"  {len(failures)} method(s) failed: {failures}")
        return 1
    print(f"  all {len(METHODS)} methods OK")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(_self_test())
