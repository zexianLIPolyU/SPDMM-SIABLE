"""
spdmm.spdmm_hybrid_rust
=======================

The "spdmm (Rust)" hybrid described in §V-A of Li, Zhang & Zhang (2025).

This module re-implements ``state_preparation`` so that the **leaf-level
unitary syntheses** (Phase 3 unitary case + Phase 4) are delegated to
qiskit's Rust-backed ``qs_decomposition`` instead of the pure-Python
``_recursive_decouple``.  The Phase-3 isometry shell (Section III-B)
is kept, but its three internal sub-syntheses also go through
``qs_decomposition``.

Cost relative to pure-Python SPDMM at n = 15:
    +7 CNOTs   (29 634 vs 29 627)
    ~2x faster compile time

Key architectural choices (matching the original v5 implementation that
produced the paper's measurements):

1.  ``_build_unitary_inverse_circuit(U, n)`` builds a sub-circuit that
    DIRECTLY implements ``U`` (no inversion needed by the caller).  We
    bit-reverse the matrix indices to bridge SPDMM's MSB-first
    convention with qiskit's LSB-first.

2.  Phase 3 (full-unitary case) and Phase 4 simply ``qc.compose(sub, ...)``
    -- a batch copy at the C-side, **not** a per-gate walk.  This is
    the critical optimisation: the slower alternative (build U†, then
    invert gate-by-gate) is what cost the previous version of this
    module ~1 s at n = 15.

3.  Phase 3 (isometry case) still uses the pure-Python isometry shell,
    but the three sub-syntheses inside it go through qsd.  We use
    ``_compose_inverse`` only for this path, since the shell itself
    emits U†.

Helpers are imported from ``spdmm.spdmm`` to avoid duplication.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import svd
from qiskit import QuantumCircuit

# Pure-Python SPDMM helpers (everything except the top-level recursion)
from .spdmm import (
    _u2_decomp, _apply_u3,
    _decouple_unitary, _transformed_zxz,
    _apply_ucrz,
    _compose_inverse,
    _Z,
)


# ---------------------------------------------------------------------------
# Bit reversal (MSB-first <-> LSB-first bridge)
# ---------------------------------------------------------------------------

def _bit_reverse_indices(n: int) -> np.ndarray:
    out = np.empty(1 << n, dtype=int)
    for k in range(1 << n):
        r = 0
        v = k
        for _ in range(n):
            r = (r << 1) | (v & 1)
            v >>= 1
        out[k] = r
    return out


# ---------------------------------------------------------------------------
# qsd leaf: build a sub-circuit that DIRECTLY implements U
# ---------------------------------------------------------------------------

def _build_unitary_inverse_circuit(U: np.ndarray, n_local: int):
    """Build a sub-circuit on ``n_local`` qubits that DIRECTLY implements
    ``U`` (no inversion needed by the caller), via qiskit's Rust-backed
    ``qs_decomposition``.

    SPDMM uses MATLAB/QCLAB big-endian convention (qubit 0 = MSB);
    qs_decomposition uses qiskit's little-endian.  We bit-reverse U's
    rows and columns to bridge: ``U_qiskit[i, j] = U[bitrev(i), bitrev(j)]``.
    """
    from qiskit.synthesis.unitary.qsd import qs_decomposition
    P = _bit_reverse_indices(n_local)
    U_qiskit = U[np.ix_(P, P)]
    return qs_decomposition(U_qiskit, opt_a1=True, opt_a2=True)


# ---------------------------------------------------------------------------
# Phase-3 isometry shell with qsd leaves
# ---------------------------------------------------------------------------

def _isometry_rust(qc: QuantumCircuit, U_dagger: np.ndarray,
                   n_local: int, base: int):
    """Synthesise the first 2^(n_local-1) columns of U_dagger using two
    multiplexors; the three sub-syntheses (VC, VB, UB at n_local - 1
    qubits) are routed to ``qs_decomposition``.  Returns Δ = I."""
    A_, B1, B2, C = _transformed_zxz(U_dagger)
    half = B1.shape[0]
    half_l = half // 2
    Zext = np.kron(_Z, np.eye(half_l))
    UC, DC, VC = _decouple_unitary(np.eye(half), C)
    UB, DB, VB = _decouple_unitary(B1 @ UC, B2 @ UC @ Zext)

    sub_VC = _build_unitary_inverse_circuit(VC, n_local - 1)
    sub_VB = _build_unitary_inverse_circuit(VB, n_local - 1)
    sub_UB = _build_unitary_inverse_circuit(UB, n_local - 1)

    qc.compose(sub_VC, qubits=range(base + 1, base + n_local), inplace=True)
    _apply_ucrz(qc, DC, n_local, "L", base)
    qc.h(base)
    qc.compose(sub_VB, qubits=range(base + 1, base + n_local), inplace=True)
    _apply_ucrz(qc, DB, n_local, "M", base)
    qc.compose(sub_UB, qubits=range(base + 1, base + n_local), inplace=True)
    qc.h(base)
    return np.eye(4, dtype=complex)


def _build_isometry_circuit_rust(U_dagger: np.ndarray, n_local: int):
    sub = QuantumCircuit(n_local)
    Delta = _isometry_rust(sub, U_dagger, n_local, base=0)
    return sub, Delta


# ---------------------------------------------------------------------------
# Top-level SPDMM recursion (Rust hybrid)
# ---------------------------------------------------------------------------

def _spdmm_recursive_rust(qc: QuantumCircuit, vector: np.ndarray, offset: int):
    """Recursive SPDMM driver with qsd leaves.  Δ = I throughout."""
    N = len(vector)
    n = int(round(np.log2(N)))

    if n == 1:
        v = vector / np.linalg.norm(vector)
        SU2 = np.column_stack([v, np.array([-v[1].conj(), v[0].conj()])])
        _, gamma, beta, delta = _u2_decomp(SU2)
        _apply_u3(qc, offset, gamma, beta, delta)
        return

    ceil_h = (n + 1) // 2
    floor_h = n // 2

    # ---- Steps 1 & 2: reshape + SVD ----
    mat = vector.reshape(1 << floor_h, 1 << ceil_h, order="F").T
    U_svd, S_vec, Vt = svd(mat, full_matrices=True)
    V_svd = Vt.conj().T
    S_mat = np.zeros((1 << ceil_h, 1 << floor_h), dtype=complex)
    for i, s in enumerate(S_vec):
        S_mat[i, i] = s

    # ---- Phase 3 sub-circuit ----
    phase3_sub_inv = None    # implements U_svd DIRECTLY
    phase3_sub_fwd = None    # implements U_svd† (needs _compose_inverse)
    phase3_uvals = None
    if ceil_h == 1:
        _, g3, b3, d3 = _u2_decomp(U_svd)
        phase3_uvals = (g3, b3, d3)
    elif ceil_h > floor_h and ceil_h > 2:
        phase3_sub_fwd, _ = _build_isometry_circuit_rust(
            U_svd.conj().T, ceil_h
        )
    else:
        phase3_sub_inv = _build_unitary_inverse_circuit(U_svd, ceil_h)

    # ---- Phase 4 sub-circuit ----
    phase4_sub_inv = None
    phase4_uvals = None
    if floor_h == 1:
        _, g4, b4, d4 = _u2_decomp(V_svd.conj())
        phase4_uvals = (g4, b4, d4)
    else:
        phase4_sub_inv = _build_unitary_inverse_circuit(V_svd.conj(), floor_h)

    # ---- Step 4: residual state for Phase 1 (no Δ migration) ----
    if n == 2:
        new_vec = np.diag(S_mat).astype(complex)
    else:
        d = min(S_mat.shape)
        new_vec = np.array([S_mat[i, i] for i in range(d)], dtype=complex)
    full = np.zeros(1 << floor_h, dtype=complex)
    full[: len(new_vec)] = new_vec

    # ---- Emit gates: Phase 1 -> 2 -> 3 -> 4 ----

    # Phase 1
    rec_offset = offset + (ceil_h - floor_h)
    _spdmm_recursive_rust(qc, full, rec_offset)

    # Phase 2: CNOT copy
    for i in range(1, floor_h + 1):
        ctrl = i + ceil_h - floor_h - 1 + offset
        targ = i + ceil_h - 1 + offset
        qc.cx(ctrl, targ)

    # Phase 3
    if phase3_uvals is not None:
        g3, b3, d3 = phase3_uvals
        _apply_u3(qc, offset, g3, b3, d3)
    elif phase3_sub_inv is not None:
        qc.compose(phase3_sub_inv,
                   qubits=range(offset, offset + ceil_h),
                   inplace=True)
    else:
        _compose_inverse(qc, phase3_sub_fwd, offset)

    # Phase 4
    if phase4_uvals is not None:
        g4, b4, d4 = phase4_uvals
        _apply_u3(qc, offset + ceil_h, g4, b4, d4)
    else:
        qc.compose(phase4_sub_inv,
                   qubits=range(offset + ceil_h, offset + n),
                   inplace=True)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def state_preparation_rust(vector) -> QuantumCircuit:
    """SPDMM with leaf-level unitary synthesis routed to qs_decomposition.

    Same interface as ``spdmm.state_preparation``; returns a
    QuantumCircuit whose output state matches ``vector`` up to a global
    phase.  Cost vs pure-Python SPDMM at n = 15: +7 CNOTs, ~2x faster.
    """
    v = np.asarray(vector, dtype=complex).reshape(-1)
    N = len(v)
    n = int(round(np.log2(N)))
    assert 1 << n == N, "len(vector) must be a power of 2"
    v = v / np.linalg.norm(v)

    # The internal recursion uses MATLAB/QCLAB MSB-first convention
    # (qubit 0 = MSB locally).  Permute the input vector to that order
    # so the output state, when read back in qiskit's LSB-first form,
    # matches the original ``vector``.
    perm = _bit_reverse_indices(n)
    v_mat = v[perm]

    qc = QuantumCircuit(n)
    _spdmm_recursive_rust(qc, v_mat, offset=0)
    return qc


class SpdmmInitializeRust:
    """qclib-style adapter for the Rust hybrid."""

    def __init__(self, state_vector, label="spdmm (Rust)"):
        v = np.asarray(state_vector, dtype=complex).reshape(-1)
        n = int(round(np.log2(len(v))))
        v = v / np.linalg.norm(v)
        self.definition = state_preparation_rust(v)
        self.num_qubits = n
        self.label = label


__all__ = ["SpdmmInitializeRust", "state_preparation_rust"]
