"""
siable.py
=========
Python port of the SIABLE (Single Ancilla Block Encoding) algorithm of
Li, Zhang, Zhang (2025), based on the MATLAB/QCLAB reference shipped
with the paper

    "Reducing C-NOT Counts for State Preparation and Block Encoding via
     Diagonal Matrix Migration".

This module exposes two public functions:

    * ``siable(matrix)``           -- single-ancilla block encoding for a
                                     general (full-rank) ``2^n x 2^n`` matrix,
                                     with the *spectral norm* as the
                                     sub-normalisation factor.  CNOT leading
                                     constant ``11/48 * 4^(n+1)``.
    * ``siable_low_rank(matrix, K)`` -- single-ancilla block encoding for a
                                     rank-``K`` matrix.  CNOT leading
                                     constant ``(K + 11/12) * 2^(n+1)``.

Both functions also return the sub-normalisation factor ``alpha = ||A||_2``,
so the user can recover ``A`` from the circuit output via
``alpha * <0_anc|U|0_anc>``.

References
----------
*   [Paper]  Li, Zhang, Zhang. *Reducing C-NOT Counts for State Preparation
             and Block Encoding via Diagonal Matrix Migration*. (2025)
*   [Krol]   Krol & Al-Ars. *Beyond Quantum Shannon Decomposition: Circuit
             Construction for n-qubit Gates Based on Block-ZXZ
             Decomposition*. Phys. Rev. Appl. 22, 034019 (2024).
*   [Vale]   Vale, Azevedo, Araujo, Araujo, da Silva. *Circuit Decomposition
             of Multicontrolled Special Unitary Single-Qubit Gates*. IEEE
             TCAD 43(3), 802-811 (2024).
*   [Iten]   Iten, Colbeck, Kukuljan, Home, Christandl. *Quantum Circuits
             for Isometries*. Phys. Rev. A 93, 032318 (2016).

Implementation overview
-----------------------
**Full-rank** (Section IV-A of the paper, Theorem 1).

Given ``A`` with ``||A||_2 = alpha``, scale ``A_norm = A / alpha`` and
SVD-decompose ``A_norm = W S V^dag``.  Define
``D = diag(exp(i * arccos(s_i)))``.  Then

    A_norm = (W D V^dag + W D^dag V^dag) / 2
           = <0| (H_anc x I) (A1 (+) A2) (H_anc x I) |0>

where ``A1 = W D V^dag`` and ``A2 = W D^dag V^dag``.  The factorisation

    diag(A1, A2) = diag(W, W) . diag(D, D^dag) . diag(V^dag, V^dag)

reduces the encoding to: (1) a Hadamard on the ancilla, (2) ``V^dag``
applied to the data, (3) a uniformly-controlled R_Z gate (the
``diag(D, D^dag)`` multiplexor), (4) ``W`` on the data, (5) a Hadamard on
the ancilla.  The two unitary syntheses use the Block-ZXZ recursive
decoupling of Krol et al.; the diagonal matrix residual ``Delta``
emerging from synthesising ``V^dag`` is migrated through the UCRZ (they
commute -- both are diagonal in the computational basis) and absorbed
into the subsequent ``W`` synthesis, saving the C-NOTs needed for one
diagonal.

**Low-rank** (Section IV-B of the paper, Theorem 4).

For a rank-``K`` matrix, only the first ``K`` columns/rows of ``W`` and
``V^dag`` are non-trivial.  We replace the full unitary synthesis by an
*isometry* synthesis (Iten 2016, column-by-column decomposition), which
costs ``K * 2^n`` C-NOTs at leading order rather than ``11/24 * 4^(n+1)``
for the full unitary.  We use qiskit's built-in ``Isometry`` class as
the isometry primitive; qiskit's recent versions (>= 0.24) use the
Vale-et-al. multicontrolled-SU(2) decomposition (``20n - 38`` CNOTs)
internally, so the resulting circuit reaches the paper's leading order.

Conventions
-----------
The reference MATLAB code uses the QCLAB convention ``qubit 0 == MSB``;
Qiskit uses ``qubit 0 == LSB``.  As in ``spdmm_siable.state_preparation``,
the algorithm is implemented internally in the MATLAB/big-endian
convention.  The public functions bit-reverse the input *matrix* (rows
and columns) on the data sub-system so the output circuit, interpreted
under Qiskit's natural ordering, block-encodes the original ``A``::

    A[i, j] == alpha * U[2i, 2j]    for the returned circuit U

i.e. the ancilla is qiskit qubit 0 and the data sub-system spans qubits
``1 .. n``.

Validation
----------
The implementation has been verified to:

1.  Match the C-NOT counts reported in Table III of the paper for n = 3..7
    (full-rank case): 9, 45, 205, 877, 3629.
2.  Block-encode arbitrary input matrices to numerical precision
    (``||A - alpha * <0|U|0>||_2 < 1e-12`` for n <= 6).
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import svd

from qiskit import QuantumCircuit
from qiskit.circuit import Gate
from qiskit.circuit.library import Isometry as _QiskitIsometry
from qiskit.quantum_info import Statevector

# Re-use the SPDMM machinery -- this is the same recursive Block-ZXZ
# decomposition used for state preparation in spdmm_siable.py.
from .spdmm import (
    _recursive_decouple,
    _apply_ucrz,
    _bit_reverse_indices,
    _transformed_zxz,
    _decouple_unitary,
    _Z,
    state_preparation,
)


# ---------------------------------------------------------------------------
# Phase-correction helper
# ---------------------------------------------------------------------------

def _fix_global_phase(qc: QuantumCircuit, A_norm: np.ndarray, n: int) -> None:
    """
    Absorb the residual global phase mismatch produced by the recursive
    Block-ZXZ synthesis into ``qc.global_phase``.

    The Block-ZXZ U(4) base case in ``spdmm_siable._recursive_decouple``
    drops a unit-modulus scalar factor at each level (the determinant
    phase ``det(U)^{1/4}`` of each 2-qubit unitary).  This is harmless
    for *state preparation* (state vectors live in projective space), but
    for block encoding it shifts every matrix entry by the same global
    phase.  The MATLAB reference (``siable.m``) tracks this phase
    analytically through ``recursive_decouple_quantum_circuit`` and
    applies ``RotationZ(0, -2*angle(global_phase))`` on the ancilla as
    the final gate.  Here we instead compute the phase a posteriori from
    a single Statevector simulation and absorb it into ``qc.global_phase``,
    which keeps the gate-count optimal.

    Parameters
    ----------
    qc : QuantumCircuit
        The (n+1)-qubit circuit produced by SIABLE so far, before phase
        correction.  Modified in place: ``qc.global_phase`` is updated.
    A_norm : (N, N) ndarray
        The alpha-normalised input matrix ``A / ||A||_2`` *in its
        original (qiskit-natural) ordering* -- NOT bit-reversed.
    n : int
        ``log2(N)``; the data sub-system size.  Total qubits = ``n + 1``.

    Cost
    ----
    One Statevector simulation of ``qc`` on a basis state.  Memory is
    ``O(2^(n+1))``; runtime is ``O(2^(n+1) * num_gates)`` plus qiskit
    overhead.  For n up to roughly 18, this is significantly cheaper
    than ``Operator(qc)`` (``O(4^(n+1))``).
    """
    N = 1 << n
    # Pick the largest-magnitude entry of A_norm for numerical robustness.
    abs_A = np.abs(A_norm)
    flat_idx = int(np.argmax(abs_A))
    i, j = flat_idx // N, flat_idx % N
    a_ij = A_norm[i, j]
    if abs(a_ij) < 1e-14:
        # A is essentially zero -- no phase to fix.
        return

    # Prepare the basis state |ancilla=0, data=j> in qiskit's natural
    # ordering (qubit 0 = ancilla = LSB; qubits 1..n = data qubits,
    # qubit 1 = LSB of j).  The flat index is 2*j.
    init = QuantumCircuit(n + 1)
    for k in range(n):
        if (j >> k) & 1:
            init.x(k + 1)
    init.compose(qc, inplace=True)

    sv = Statevector.from_instruction(init).data
    # The block-encoding row index is 2*i (ancilla=0 component, data=i).
    u_2i_2j = complex(sv[2 * i])

    # We want:  u_2i_2j == A_norm[i, j].
    # Currently we have:  u_2i_2j == e^{i*phi} * A_norm[i, j].
    # Setting qc.global_phase -= phi multiplies every matrix entry
    # by e^{-i*phi}, exactly cancelling the spurious phase.
    phase = u_2i_2j / a_ij
    if abs(phase) < 1e-14:
        return
    phi = float(np.angle(phase))
    qc.global_phase = float(qc.global_phase) - phi


# ---------------------------------------------------------------------------
# Bit-reversal helpers for matrices
# ---------------------------------------------------------------------------

def _bit_reverse_matrix(A: np.ndarray, n: int) -> np.ndarray:
    """
    Bit-reverse the rows AND columns of an ``n``-qubit operator ``A``
    (size ``2^n x 2^n``).  This converts between Qiskit's natural
    little-endian ordering and the MATLAB/QCLAB big-endian ordering used
    internally by the SPDMM-style code.

    Equivalent to ``P @ A @ P`` where ``P`` is the bit-reversal
    permutation matrix on n qubits (``P[i, j] = delta(j, bitrev_n(i))``).
    """
    P = _bit_reverse_indices(n)
    return A[np.ix_(P, P)]


# ---------------------------------------------------------------------------
# Explicit C-NOT cost functions
# ---------------------------------------------------------------------------
#
# These are exact formulas (or empirical look-ups) for the various
# building blocks used by SIABLE.  They let us *predict* the C-NOT count
# of an isometry *before* synthesising the circuit, so the low-rank
# driver can pick the cheapest method without paying the (sometimes
# expensive) build cost.
#
# Conventions used below:
#   * ``m`` denotes the number of qubits the routine acts on.
#   * ``K`` denotes the number of columns of an isometry / a power of 2.
#
# Derivations
# -----------
# 1.  Uniformly-controlled R_Z (``_apply_ucrz``) on ``m`` qubits:
#
#         UCRZ_L(m) = UCRZ_R(m) = 2^{m-1} - 1
#         UCRZ_M(m) = 2^{m-1}                          (exact)
#
# 2.  Block-ZXZ recursive unitary synthesis ``_recursive_decouple`` on
#     ``m`` qubits (``is_last_U=False`` -- i.e. the 2-CNOT KAK base case):
#
#         f(2) = 2,
#         f(m) = 4 f(m-1) + 3 * 2^{m-1} - 2     (m >= 3)
#
#     Closed form (verified against build counts up to m=8):
#
#         f(m) = (11/24) * 4^m - 3 * 2^{m-1} + 2/3
#
# 3.  SPDMM half-isometry ``_isometry`` on ``m`` qubits (gives the first
#     ``2^{m-1}`` columns of an ``m``-qubit unitary, valid for m >= 3):
#
#         h(m) = 3 * f(m-1) + UCRZ_L(m) + UCRZ_M(m)
#              = 3 * f(m-1) + 2^m - 1
#              = (11/32) * 4^m - (5/4) * 2^m + 1
#
# 4.  SPDMM state preparation (``state_preparation``) on ``m`` qubits
#     (i.e. a 1-column isometry from |0> to a 2^m-dim state):
#
#     Empirical values measured for m = 2..10:
#
#         m: 2  3  4  5   6   7   8    9    10
#         #: 1  3  7  18  42  93  199  411  837
#
#     We store these in a lookup table; the few-millisecond build for any
#     m not in the table is also acceptable (fallback below).
#
# 5.  qiskit ``Isometry`` on ``m`` qubits with ``K`` (=power of 2)
#     columns: leading-order ``~ K * 2^m``.  We don't have a simple
#     closed-form across all (m, K), so we BUILD a probe and count
#     C-NOTs.  This is fast (qiskit's Isometry synthesis itself is
#     cheap; only its execution / Statevector simulation is expensive).
#

# Look-up table for SPDMM state preparation C-NOT counts, populated at
# import time the first time the cost function is called.  Keys are the
# number of qubits ``m``.
_SPDMM_STATE_PREP_CNOTS = {}


def _cnot_cost_recursive_decouple(m: int) -> int:
    """
    Exact C-NOT count for an ``m``-qubit Block-ZXZ unitary synthesis
    (``_recursive_decouple`` with ``is_last_U=False``).

    Returns 0 for ``m <= 1``; the closed-form formula for ``m >= 2`` is

        f(m) = (11/24) * 4^m - 3 * 2^{m-1} + 2/3.

    The expression is always an integer for integer m >= 2.
    """
    if m <= 1:
        return 0
    val = (11 / 24) * (4 ** m) - 3 * (2 ** (m - 1)) + 2 / 3
    return int(round(val))


def _cnot_cost_spdmm_half_iso(m: int) -> int:
    """
    Exact C-NOT count for an ``m``-qubit SPDMM half-isometry as used by
    :func:`_siable_low_rank_spdmm_half` -- i.e. with ``is_last_U=True``
    at the deepest base case, so the residual Delta is identity and no
    diagonal needs to be migrated through the UCRZ.

    Valid for ``m >= 3``; returns ``+inf`` otherwise to mark "not
    applicable" (the half-isometry routine internally recurses to
    ``m - 1`` qubits which requires ``m - 1 >= 2``).

    Closed form

        h_T(m) = (11/32) * 4^m - (5/4) * 2^m + 2,

    which is ``h(m) + 1`` where ``h(m)`` is the standard half-iso
    cost (with ``is_last_U=False``).  The extra ``+1`` comes from the
    3-CNOT KAK base case (vs the 2-CNOT KAK used when leakage is
    allowed).
    """
    if m < 3:
        return int(10 ** 18)  # sentinel "not applicable"
    val = (11 / 32) * (4 ** m) - (5 / 4) * (2 ** m) + 2
    return int(round(val))


def _cnot_cost_spdmm_state_prep(m: int) -> int:
    """
    C-NOT count for ``spdmm_siable.state_preparation`` on ``m`` qubits.

    Uses a cached table for the common range ``m = 1..12``; if a value
    is missing, builds a single probe state and counts ``cx`` gates.
    Random vectors and structured ones produce the same count, so this
    is a property of ``m`` alone.
    """
    if m <= 0:
        return 0
    if m == 1:
        return 0
    if m in _SPDMM_STATE_PREP_CNOTS:
        return _SPDMM_STATE_PREP_CNOTS[m]
    # Build a small probe and count.
    rng = np.random.default_rng(0)
    v = rng.standard_normal(1 << m) + 1j * rng.standard_normal(1 << m)
    v = v / np.linalg.norm(v)
    qc_probe = state_preparation(v)
    cnots = qc_probe.count_ops().get("cx", 0)
    _SPDMM_STATE_PREP_CNOTS[m] = int(cnots)
    return int(cnots)


def _cnot_cost_qiskit_isometry(m: int, K: int) -> int:
    """
    C-NOT count for qiskit's ``Isometry`` synthesis with ``K`` columns
    on ``m`` qubits (``K`` must be a power of 2, ``K <= 2^m``).

    qiskit chooses its multi-controlled-SU(2) implementation based on
    the random matrix structure only mildly (the choice of axis-aligned
    or Iten-style decomposition depends on the column shape but not its
    random entries).  We therefore build a probe isometry of the right
    shape and count C-NOTs; the result is essentially independent of
    the probe and is cached on ``(m, K)``.
    """
    if K == 1:
        # Closed form for state preparation via qiskit Isometry.
        return (1 << m) - m - 1
    key = ("qiskit_iso", m, K)
    if key in _CNOT_COST_CACHE:
        return _CNOT_COST_CACHE[key]
    rng = np.random.default_rng(0)
    N = 1 << m
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    Q, _ = np.linalg.qr(A)
    V_iso = Q[:, :K]
    iso_gate = _QiskitIsometry(V_iso, 0, 0)
    qc_probe = QuantumCircuit(m)
    qc_probe.append(iso_gate, list(range(m)))
    qc_probe = qc_probe.decompose(reps=6)
    cnots = int(qc_probe.count_ops().get("cx", 0))
    _CNOT_COST_CACHE[key] = cnots
    return cnots


def _cnot_cost_full_rank(m: int) -> int:
    """
    Exact total C-NOT count for the full-rank ``siable()`` construction
    on an ``(m+1)``-qubit circuit with ``m`` data qubits:

        full(m) = 2 * f(m) + 1 + 2^m

    where ``f(m) = (11/24)*4^m - 3*2^(m-1) + 2/3`` is the Block-ZXZ
    unitary-synthesis cost.  The ``+1`` is the 3-CNOT KAK base case
    used by the W-side synthesis (``is_last_U=True``).  The ``2^m``
    term is ``UCRZ_M(m+1)``.

    Matches paper Table III exactly: full(2)=9, full(3)=45, full(4)=205,
    full(5)=877, full(6)=3629.
    """
    if m < 2:
        return 0
    f_m = _cnot_cost_recursive_decouple(m)
    return 2 * f_m + 1 + (1 << m)


# Generic CNOT-cost cache for non-trivial builds.
#
# Pre-populated with known qiskit ``Isometry`` C-NOT counts for the
# (m_qubits, K_columns) pairs commonly encountered by ``siable_low_rank``.
# Without these, the first call for each new (m, K) pair has to BUILD a
# probe qiskit Isometry on ``m`` qubits and recursively decompose it
# (``decompose(reps=6)``), which is ``O(K * 4^m)`` time -- the
# ``(m=8, K=16)`` probe alone takes ~2 s, and the ``(m=7, K=128)``
# probe ~10 s.  By baking the answer in we avoid that build entirely.
#
# Values verified against qiskit 2.4.1 on Linux/x86_64; they are
# implementation-defined but stable across qiskit minor versions.  Cache
# misses (e.g. for very large m, or unusual K) fall back to the on-the-
# fly probe build below.
_CNOT_COST_CACHE = {
    ("qiskit_iso", 2, 2): 3,
    ("qiskit_iso", 2, 4): 6,
    ("qiskit_iso", 3, 2): 10,
    ("qiskit_iso", 3, 4): 24,
    ("qiskit_iso", 3, 8): 41,
    ("qiskit_iso", 4, 2): 25,
    ("qiskit_iso", 4, 4): 57,
    ("qiskit_iso", 4, 8): 122,
    ("qiskit_iso", 4, 16): 218,
    ("qiskit_iso", 5, 2): 56,
    ("qiskit_iso", 5, 4): 122,
    ("qiskit_iso", 5, 8): 261,
    ("qiskit_iso", 5, 16): 546,
    ("qiskit_iso", 5, 32): 1025,
    ("qiskit_iso", 6, 2): 119,
    ("qiskit_iso", 6, 4): 251,
    ("qiskit_iso", 6, 8): 528,
    ("qiskit_iso", 6, 16): 1107,
    ("qiskit_iso", 6, 32): 2300,
    ("qiskit_iso", 6, 64): 4474,
    ("qiskit_iso", 7, 2): 246,
    ("qiskit_iso", 7, 4): 508,
    ("qiskit_iso", 7, 8): 1051,
    ("qiskit_iso", 7, 16): 2180,
    ("qiskit_iso", 7, 32): 4527,
    ("qiskit_iso", 7, 64): 9372,
    ("qiskit_iso", 7, 128): 18653,
    ("qiskit_iso", 8, 2): 501,
    ("qiskit_iso", 8, 4): 1021,
    ("qiskit_iso", 8, 8): 2086,
    ("qiskit_iso", 8, 16): 4277,
}


# ---------------------------------------------------------------------------
# SPDMM half-isometry with controllable final ``is_last_U`` flag
# ---------------------------------------------------------------------------

def _isometry_with_last_U(qc: QuantumCircuit, U_dagger: np.ndarray,
                          n_local: int, base: int, Delta: np.ndarray,
                          is_last_U: bool):
    """
    Same as ``spdmm_siable._isometry`` -- builds a half-isometry sub-
    circuit on ``n_local`` qubits whose first ``2^{n_local-1}`` columns
    realise the first ``2^{n_local-1}`` columns of ``U_dagger`` (up to a
    4x4 residual diagonal on the bottom two qubits) -- BUT with an
    extra ``is_last_U`` flag controlling whether the final
    ``_recursive_decouple`` call uses the 2-CNOT KAK (``False``, leaks
    a diagonal) or the 3-CNOT KAK (``True``, clean exit).

    Parameters
    ----------
    qc : QuantumCircuit
    U_dagger : (2^n_local, 2^n_local) complex ndarray
        The "target" unitary; only its first ``2^{n_local-1}`` columns
        are pinned down by the synthesis.
    n_local : int
        ``>= 3``.  (For ``n_local == 2`` the half-isometry has only one
        column and is a single-state preparation; we never need that
        path in this file.)
    base : int
        Qubit offset where this sub-circuit begins.
    Delta : (4, 4) complex ndarray
        Input residual diagonal coming from a preceding synthesis
        (identity if this is the first synthesis).
    is_last_U : bool
        If True, the third recursive call uses ``is_last_U=True`` and
        no diagonal leaks at the bottom; if False, a 4x4 diagonal Delta
        leaks (returned).  Mirrors the convention used by ``siable``.

    Returns
    -------
    Delta_out : (4, 4) complex ndarray
        Identity (up to numerical noise) when ``is_last_U=True``, else
        the leaked diagonal.
    """
    A_, B1, B2, C = _transformed_zxz(U_dagger)
    half = B1.shape[0]
    half_l = half // 2
    Zext = np.kron(_Z, np.eye(half_l))
    UC, DC, VC = _decouple_unitary(np.eye(half), C)
    UB, DB, VB = _decouple_unitary(B1 @ UC, B2 @ UC @ Zext)

    Delta = _recursive_decouple(qc, VC, n_local - 1, base + 1, Delta, False)
    _apply_ucrz(qc, DC, n_local, "L", base)
    qc.h(base)
    Delta = _recursive_decouple(qc, VB, n_local - 1, base + 1, Delta, False)
    _apply_ucrz(qc, DB, n_local, "M", base)
    Delta = _recursive_decouple(qc, UB, n_local - 1, base + 1, Delta, is_last_U)
    qc.h(base)
    return Delta


# ---------------------------------------------------------------------------
# Public API: full-rank single-ancilla block encoding
# ---------------------------------------------------------------------------

def siable(
    matrix: np.ndarray,
    *,
    return_info: bool = False,
):
    """
    Build a single-ancilla block encoding of ``matrix`` using the SIABLE
    construction of the paper (Theorem 1).

    Parameters
    ----------
    matrix : (N, N) array_like
        The matrix ``A`` to block-encode.  ``N`` must be a power of 2.
        Can be complex; need not be unitary.
    return_info : bool, optional
        If True, also return a ``dict`` containing CNOT count and other
        bookkeeping.  Default False.

    Returns
    -------
    qc : qiskit.QuantumCircuit
        A circuit on ``log2(N) + 1`` qubits.  Qubit 0 is the ancilla;
        qubits ``1 .. log2(N)`` hold the data sub-system.  The block
        encoding relation is::

            A[i, j]  ==  alpha * Statevector @ qc @ (...)
                     ==  alpha * Operator(qc).data[2*i, 2*j]

        for all ``i, j`` in ``0 .. N - 1`` (up to numerical error).
    alpha : float
        Sub-normalisation factor (= spectral norm ``||A||_2``).
    info : dict, optional
        Returned only if ``return_info=True``.  Currently has:

        ``nCNOT``  --  exact CNOT count of the returned ``qc``.

    Notes
    -----
    Implementation closely follows ``siable.m`` in the paper's reference
    package:

        1.  Subnormalise:  ``A_norm = A / ||A||_2``.
        2.  SVD:           ``A_norm = W @ Sigma @ V^dag``.
        3.  Build D:       ``D[i, i] = exp(i * arccos(Sigma[i, i]))``.
        4.  Emit:
              H(0)
              V_dag on qubits 1..n              [up to diagonal Delta_V]
              UCRZ(D) with target=0, ctrls=1..n [commutes through Delta_V]
              W . Delta_V^dag on qubits 1..n
              H(0)
              global-phase Rz on qubit 0
    """
    A = np.asarray(matrix, dtype=complex)
    N, N2 = A.shape
    if N != N2:
        raise ValueError("matrix must be square")
    n = int(round(np.log2(N)))
    if (1 << n) != N:
        raise ValueError("matrix dimension must be a power of 2")
    if n < 2:
        raise ValueError("siable requires n >= 2 (a >= 4x4 matrix)")

    # ----- Step 1: subnormalisation -----
    alpha = float(np.linalg.norm(A, ord=2))
    if alpha == 0.0:
        # Edge case: zero matrix.  Return identity-like circuit; alpha=1
        # is the natural convention.
        qc = QuantumCircuit(n + 1)
        if return_info:
            return qc, 1.0, {"nCNOT": 0}
        return qc, 1.0

    A_norm = A / alpha

    # Convert to MATLAB/QCLAB endianness on the data qubits (n qubits).
    # The ancilla is qubit 0 of the (n+1)-qubit circuit in BOTH the
    # qiskit and the internal-MATLAB pictures -- bit-reversal only acts
    # on the data sub-system.
    A_int = _bit_reverse_matrix(A_norm, n)

    # ----- Step 2: SVD of the (sub-normalised, bit-reversed) matrix -----
    W, S_vec, Vh = svd(A_int)
    # numerical safety: clamp singular values to [0, 1] before arccos
    s_clamped = np.clip(S_vec, 0.0, 1.0)
    D_diag = np.exp(1j * np.arccos(s_clamped))  # 1D array, length 2^n
    D = np.diag(D_diag)

    # ----- Step 3: build the (n+1)-qubit circuit -----
    qc = QuantumCircuit(n + 1)
    # H on ancilla
    qc.h(0)

    # V^dag on the n data qubits (base=1 in MATLAB convention).
    # is_last_U=False so that the residual 4x4 diagonal Delta leaks out
    # and we can absorb it into the subsequent W synthesis.
    V_dag = Vh  # already the conjugate-transpose
    Delta = _recursive_decouple(
        qc, V_dag, n_local=n, base=1,
        Delta=np.eye(4, dtype=complex), is_last_U=False,
    )

    # UCRZ(D) -- target = qubit 0 (ancilla), controls = qubits 1..n.
    # This diagonal multiplexor commutes with Delta (both diagonal in
    # the computational basis), so Delta is unaffected by inserting
    # UCRZ between the V^dag synthesis and the W synthesis.
    _apply_ucrz(qc, D, n_local=n + 1, mode="M", base=0)

    # W absorbs Delta^dag (the leaked diagonal from V^dag) into its
    # synthesis.  Pass Delta to _recursive_decouple as the "initial"
    # residual -- the recursion threads it through to the first U(4)
    # base case, which absorbs it.  is_last_U=True so no Delta leaks
    # out (we want a clean encoding here).
    _recursive_decouple(
        qc, W, n_local=n, base=1,
        Delta=Delta, is_last_U=True,
    )

    # H on ancilla
    qc.h(0)

    # Phase correction: the recursive Block-ZXZ synthesis drops a
    # global phase per U(4) base case.  Absorb the residual phase into
    # qc.global_phase so the block encoding has the correct sign.
    _fix_global_phase(qc, A_norm, n)

    if return_info:
        info = {"nCNOT": qc.count_ops().get("cx", 0)}
        return qc, alpha, info
    return qc, alpha


# ---------------------------------------------------------------------------
# Helper: isometry synthesis from 2^k columns of a 2^n x 2^n matrix
# ---------------------------------------------------------------------------

def _isometry_first_k_columns(
    V_full: np.ndarray, K: int, n: int,
) -> QuantumCircuit:
    """
    Build a qiskit circuit on ``n`` qubits whose action on
    ``|0>, |1>, ..., |K-1>`` (qiskit-natural ordering) produces the
    first ``K`` columns of ``V_full``.

    ``V_full`` must be expressed in qiskit-natural ordering -- i.e.
    ``V_full[:, j]`` is the desired output state when the input is
    ``|j>_qiskit``.

    Internally:
      *   Take the first ``K`` columns ``V_iso = V_full[:, :K]`` -- this
          is a ``2^n x K`` isometry (orthonormal columns).
      *   Hand to qiskit's ``Isometry`` class, which implements the
          column-by-column decomposition of Iten et al., Phys. Rev. A
          93, 032318 (2016).  Modern qiskit (>= 0.24.1) uses the
          Vale-et-al. multicontrolled-SU(2) optimisation (Phys. Rev. A
          --, 20n - 38 CNOTs per multi-controlled gate) under the hood.

    Returns
    -------
    sub : QuantumCircuit
        An ``n``-qubit circuit ``C`` such that, in qiskit's natural
        ordering, ``C |j> = V_full[:, j]`` for ``j = 0, 1, ..., K-1``.
    """
    V_iso = V_full[:, :K]
    iso_gate = _QiskitIsometry(V_iso, 0, 0)
    sub = QuantumCircuit(n)
    sub.append(iso_gate, list(range(n)))
    # Pre-decompose so we get a circuit at the (u, cx) level.  Without
    # this, ``sub`` is a single opaque Isometry gate that the qiskit
    # transpiler would otherwise unfold lazily.  Pre-decomposing here
    # lets the rest of the SIABLE assembly inspect the C-NOT count via
    # ``qc.count_ops()``.
    sub = sub.decompose(reps=6)
    return sub


# ---------------------------------------------------------------------------
# Public API: rank-K single-ancilla block encoding
# ---------------------------------------------------------------------------

def siable_low_rank(
    matrix: np.ndarray,
    K: int,
    *,
    method: str = "auto",
    return_info: bool = False,
):
    """
    Build a single-ancilla block encoding of a rank-``K`` matrix
    ``matrix`` using the low-rank SIABLE construction (Theorem 4 of the
    paper), optimised to pick the cheapest construction automatically.

    For a matrix ``A`` of shape ``2^n x 2^n`` and target rank ``K``
    (``1 <= K <= 2^n``), the rank-``K`` truncation of ``A`` is

        A_K = W @ Sigma_K @ V^dag

    where ``Sigma_K`` is the top-``K`` truncation of the singular-value
    matrix.  The block encoding subnormalises by ``||A||_2``.

    Parameters
    ----------
    matrix : (N, N) array_like
    K : int
        Rank to encode.  Must satisfy ``K >= 1``.  Values of ``K`` larger
        than ``N`` are silently clamped to ``N`` (since the rank-``K``
        truncation of an ``N x N`` matrix is just ``A`` itself once
        ``K >= N``).
    method : {'auto', 'state_prep', 'spdmm_half', 'isometry', 'full_rank'}
        Which construction to use.

        - ``'auto'`` (default): explicit C-NOT-cost estimates pick the
          cheapest among **four** candidates -- ``'state_prep'``
          (K=1 only), ``'spdmm_half'`` (K_pad <= 2^(n-1), n >= 3),
          ``'isometry'``, and ``'full_rank'``.

          **Note**: when ``'full_rank'`` wins the comparison, the
          returned circuit block-encodes ``A`` directly, NOT the
          rank-``K`` truncation ``A_K``.  Since ``A`` is a strictly
          better approximation of itself than ``A_K`` is, this is
          beneficial whenever cost permits.  If the caller specifically
          requires the rank-``K`` truncation regardless of cost, pass
          ``method='isometry'`` (or another non-``'full_rank'``
          option).

        - ``'state_prep'``: force the K=1 fast path; raises if ``K > 1``.
        - ``'spdmm_half'``: force the SPDMM half-isometry path; raises
          if ``K_pad > 2^(n-1)`` or ``n < 3``.
        - ``'isometry'``: force qiskit's ``Isometry``-based path.
        - ``'full_rank'``: force the full-rank ``siable()`` construction
          (encodes ``A``, not ``A_K``).
    return_info : bool, optional

    Returns
    -------
    qc : qiskit.QuantumCircuit
        A circuit on ``n + 1`` qubits encoding either ``A_K / alpha`` or
        ``A / alpha`` (see notes on ``method='auto'`` above).
    alpha : float
        ``||A||_2``.
    info : dict, optional
        Includes ``method`` (chosen path), ``predicted_cost`` (the full
        cost table considered), ``nCNOT`` (actual C-NOT count of the
        built circuit), and ``encodes`` (either ``'A_K'`` or ``'A'``).

    Notes
    -----
    Total leading-term C-NOT count of the qiskit-Isometry path:
    ``(K + 11/12) * 2^(n+1)``  (Theorem 4 of the paper).  This grows
    linearly in ``K``; at ``K ~ N/4`` it reaches the full-rank cost
    ``~(11/12) * 4^(n+1)``, at which point ``'full_rank'`` wins.
    """
    A = np.asarray(matrix, dtype=complex)
    N, N2 = A.shape
    if N != N2:
        raise ValueError("matrix must be square")
    n = int(round(np.log2(N)))
    if (1 << n) != N:
        raise ValueError("matrix dimension must be a power of 2")
    if n < 2:
        raise ValueError("siable_low_rank requires n >= 2 (a >= 4x4 matrix)")
    if K < 1:
        raise ValueError(f"K must satisfy K >= 1, got K={K}")
    # K > N has the same meaning as K = N: the rank-K truncation of an
    # N x N matrix is just A itself once K reaches N.  Silently clamp.
    K = min(K, N)

    valid_methods = {"auto", "state_prep", "spdmm_half", "isometry", "full_rank"}
    if method not in valid_methods:
        raise ValueError(
            f"method must be one of {sorted(valid_methods)}, got {method!r}"
        )

    K_pad = 1 << int(np.ceil(np.log2(max(K, 1))))
    K_pad = max(1, min(K_pad, N))

    # ----- method auto-dispatch via explicit total-C-NOT cost --------------
    # All entries below are TOTAL C-NOT counts (including the shared 2^n
    # UCRZ_M overhead for the low-rank paths) so they can be compared
    # head-to-head with the full-rank baseline.
    cost_table = {}
    ucrz_m = 1 << n  # UCRZ_M(n+1) -- shared by all low-rank paths

    # Full-rank baseline -- a viable option at any K (encodes A, not A_K).
    cost_table["full_rank"] = _cnot_cost_full_rank(n)

    if K == 1:
        cost_table["state_prep"] = 2 * _cnot_cost_spdmm_state_prep(n) + ucrz_m
        cost_table["isometry"] = 2 * _cnot_cost_qiskit_isometry(n, 1) + ucrz_m
    else:
        cost_table["isometry"] = (
            2 * _cnot_cost_qiskit_isometry(n, K_pad) + ucrz_m
        )
        if K_pad <= (1 << (n - 1)) and n >= 3:
            cost_table["spdmm_half"] = (
                2 * _cnot_cost_spdmm_half_iso(n) + ucrz_m
            )

    if method == "auto":
        method = min(cost_table, key=cost_table.get)
    else:
        # Validate the user's manual choice.
        if method == "state_prep" and K != 1:
            raise ValueError("method='state_prep' is only valid for K == 1")
        if method == "spdmm_half":
            if n < 3:
                raise ValueError("method='spdmm_half' requires n >= 3")
            if K_pad > (1 << (n - 1)):
                raise ValueError(
                    "method='spdmm_half' requires K_pad <= 2^(n-1); "
                    f"got K_pad={K_pad}, 2^(n-1)={1<<(n-1)}"
                )

    # ----- dispatch ------------------------------------------------------
    if method == "state_prep":
        qc, alpha = _siable_low_rank_state_prep(A)
        encodes = "A_K"
    elif method == "spdmm_half":
        qc, alpha = _siable_low_rank_spdmm_half(A, K)
        encodes = "A_K"
    elif method == "full_rank":
        # Use the full-rank construction directly.  This block-encodes
        # A itself (a strictly better approximation than A_K).
        qc, alpha = siable(A)
        encodes = "A"
    else:  # 'isometry'
        qc, alpha = _siable_low_rank_isometry(A, K)
        encodes = "A_K"

    if return_info:
        info = {
            "nCNOT": qc.count_ops().get("cx", 0),
            "rank": K,
            "method": method,
            "predicted_cost": cost_table,
            "encodes": encodes,
        }
        return qc, alpha, info
    return qc, alpha


# ---------------------------------------------------------------------------
# Low-rank SIABLE implementations (three paths)
# ---------------------------------------------------------------------------

def _siable_low_rank_state_prep(A: np.ndarray):
    """
    K = 1 fast path: replace qiskit ``Isometry`` with
    ``spdmm_siable.state_preparation`` for both V- and W-syntheses.

    Stays in qiskit ordering throughout (state_preparation returns a
    qiskit-ordered circuit).  Gives the *exact* C-NOT counts of the
    paper's Table V for ``n_total = 3..7`` and beats it for ``n = 8``.
    """
    N = A.shape[0]
    n = int(round(np.log2(N)))

    alpha = float(np.linalg.norm(A, ord=2))
    if alpha == 0.0:
        return QuantumCircuit(n + 1), 1.0

    A_norm = A / alpha

    # SVD in qiskit ordering (no bit-reversal in this path).
    W_q, S_vec, Vh_q = svd(A_norm)
    V_q = Vh_q.conj().T
    s_clamped = np.clip(S_vec, 0.0, 1.0)

    # Diagonal phases: only the K=1 entry is "live"; the rest are i so
    # that (D + D†)/2 = 0 on those columns.
    D_qiskit = np.empty(N, dtype=complex)
    D_qiskit[0] = np.exp(1j * np.arccos(s_clamped[0]))
    D_qiskit[1:] = 1j

    # Build the two state-prep sub-circuits.  ``state_preparation(v)``
    # returns a qiskit-ordered circuit ``C`` with ``C |0> = v`` (up to a
    # global phase that ``_fix_global_phase`` reconciles at the end).
    C_V = state_preparation(V_q[:, 0])
    C_W = state_preparation(W_q[:, 0])

    qc = QuantumCircuit(n + 1)

    # H on ancilla.
    qc.h(0)

    # Apply V^dag to data:  C_V.inverse() maps V[:,0] -> |0>.
    qc.compose(C_V.inverse(), qubits=list(range(1, n + 1)), inplace=True)

    # Uniformly-controlled R_Z(D) on (anc, data).  ``_apply_ucrz`` uses
    # MATLAB convention (data qubit 1 == MSB) but our diagonal is in
    # qiskit ordering; bit-reverse the 1D array to bridge the gap.
    P = _bit_reverse_indices(n)
    D_mat = np.diag(D_qiskit[P])
    _apply_ucrz(qc, D_mat, n_local=n + 1, mode="M", base=0)

    # Apply W to data:  C_W maps |0> -> W[:,0].
    qc.compose(C_W, qubits=list(range(1, n + 1)), inplace=True)

    # H on ancilla.
    qc.h(0)

    # Phase correction against the rank-1 truncation.
    Sigma_1 = np.zeros(N)
    Sigma_1[0] = s_clamped[0]
    A_1_qiskit = W_q @ np.diag(Sigma_1) @ Vh_q
    _fix_global_phase(qc, A_1_qiskit, n)

    return qc, alpha


def _siable_low_rank_spdmm_half(A: np.ndarray, K: int):
    """
    Use the SPDMM half-isometry (``spdmm_siable._isometry``) for both
    V- and W-syntheses.  Valid when ``K_pad <= 2^(n-1)`` and ``n >= 3``;
    pads ``K`` up to the half-isometry's fixed column count
    ``K_eff = 2^{n-1}`` by setting the corresponding singular values to
    zero (D[k] = i in those slots).

    Implementation
    --------------
    The SPDMM ``_isometry`` synthesises a sub-circuit whose MATLAB-
    ordered matrix has the property:

        M_mat[:N/2, :]  =  c * (Delta on bot 2 qubits) * U_dagger[:N/2, :]

    i.e. its **first 2^{n-1} rows** match ``U_dagger``'s first 2^{n-1}
    rows up to a global phase ``c`` and a 4x4 residual diagonal
    ``Delta`` on the two least-significant (in MATLAB convention) data
    qubits.  Setting ``is_last_U=True`` collapses ``Delta = I`` at the
    cost of one extra CNOT (3-CNOT KAK at the deepest base case),
    giving exactly ``M_mat[:N/2, :] = c * U_dagger[:N/2, :]``.

    Putting this together with the SIABLE protocol amplitude
    derivation:

        amplitude[i, j]  =  sum_{k<K} sigma_k * M_W[i, k] * M_V[k, j]

    we need ``M_V``'s first ``K`` rows to match ``V_m^dag``'s first
    ``K`` rows, and ``M_W``'s first ``K`` columns to match ``W_m``'s
    first ``K`` columns (both up to a global phase, handled by
    ``_fix_global_phase``).  The natural way to get this from the
    SPDMM half-isometry is:

      * pass ``V_m^dag`` to ``_isometry_with_last_U`` for sub_V, then
        apply sub_V **directly** to the data; its first 2^{n-1} rows
        match ``V_m^dag``'s first 2^{n-1} rows;
      * pass ``W_m^dag`` to ``_isometry_with_last_U`` for sub_W, then
        apply sub_W's **inverse** to the data; the inverse's first
        2^{n-1} columns match ``W_m``'s first 2^{n-1} columns.

    No Delta migration through the UCRZ is needed because both Deltas
    are identity (``is_last_U=True``).

    Cost
    ----
    Each ``_isometry_with_last_U`` with ``is_last_U=True`` costs
    ``h(n) + 1`` C-NOTs where ``h(n) = (11/32)*4^n - (5/4)*2^n + 1``
    (the extra ``+1`` is the 3-CNOT KAK base case).  Plus
    ``UCRZ_M(n+1) = 2^n`` for the middle UCRZ, total

        2*(h(n) + 1) + 2^n  =  (11/16)*4^n - (3/2)*2^n + 4.

    For ``n >= 4`` and ``K`` near ``2^{n-1}`` this is roughly half the
    qiskit-Isometry path's count.
    """
    N = A.shape[0]
    n = int(round(np.log2(N)))

    alpha = float(np.linalg.norm(A, ord=2))
    if alpha == 0.0:
        return QuantumCircuit(n + 1), 1.0

    A_norm = A / alpha
    # Convert to MATLAB ordering on the data sub-system.
    A_int = _bit_reverse_matrix(A_norm, n)

    # SVD in MATLAB ordering.
    W_m, S_vec, Vh_m = svd(A_int)
    s_clamped = np.clip(S_vec, 0.0, 1.0)

    # Diagonal phases, MATLAB-ordered (length 2^n).
    # D[k] = exp(i*arccos(sigma_k)) for k < K (genuine singular values),
    # D[k] = i otherwise (zero singular value -- kills the contribution
    # via (D[k] + D[k]^*) / 2 = cos(arccos(0)) = 0).
    D_diag = np.empty(N, dtype=complex)
    for k in range(K):
        D_diag[k] = np.exp(1j * np.arccos(s_clamped[k]))
    for k in range(K, N):
        D_diag[k] = 1j
    D_mat = np.diag(D_diag)

    # Build the two half-isometry sub-circuits independently (each is
    # a standalone n-qubit circuit).  ``is_last_U=True`` so the
    # residual Delta is identity in both -- no migration needed.
    sub_V = QuantumCircuit(n)
    _isometry_with_last_U(
        sub_V, Vh_m, n_local=n, base=0,
        Delta=np.eye(4, dtype=complex), is_last_U=True,
    )

    sub_W = QuantumCircuit(n)
    _isometry_with_last_U(
        sub_W, W_m.conj().T, n_local=n, base=0,
        Delta=np.eye(4, dtype=complex), is_last_U=True,
    )

    # Assemble the (n+1)-qubit SIABLE circuit.
    qc = QuantumCircuit(n + 1)

    # H on ancilla.
    qc.h(0)

    # Apply sub_V directly to data qubits 1..n.  Its first 2^{n-1}
    # MATLAB rows match V_m^dag's first 2^{n-1} rows up to a global
    # phase (the K-relevant rows, since K <= 2^{n-1}).
    qc.compose(sub_V, qubits=list(range(1, n + 1)), inplace=True)

    # UCRZ(diag(D, D^dag)) on (anc, data) -- target = qubit 0 = ancilla,
    # controls = data qubits 1..n.  In MATLAB convention (where qubit 1
    # is the data MSB).
    _apply_ucrz(qc, D_mat, n_local=n + 1, mode="M", base=0)

    # Apply sub_W^dag (inverse) to data qubits 1..n.  This gives a
    # matrix whose first 2^{n-1} columns match W_m's first 2^{n-1}
    # columns (up to global phase).
    qc.compose(sub_W.inverse(), qubits=list(range(1, n + 1)), inplace=True)

    # H on ancilla.
    qc.h(0)

    # Phase correction against the rank-K truncation in qiskit ordering.
    Sigma_K = np.zeros(N)
    Sigma_K[:K] = s_clamped[:K]
    # Rank-K truncation in MATLAB ordering = W_m @ diag(Sigma_K) @ Vh_m.
    A_int_K = W_m @ np.diag(Sigma_K) @ Vh_m
    A_K_qiskit = _bit_reverse_matrix(A_int_K, n)
    _fix_global_phase(qc, A_K_qiskit, n)

    return qc, alpha


def _siable_low_rank_isometry(A: np.ndarray, K: int):
    """
    Original qiskit-``Isometry`` based path (no bit-reversal of the
    matrix; everything in qiskit-natural ordering except the diagonal
    array fed to ``_apply_ucrz``).
    """
    N = A.shape[0]
    n = int(round(np.log2(N)))

    alpha = float(np.linalg.norm(A, ord=2))
    if alpha == 0.0:
        return QuantumCircuit(n + 1), 1.0

    A_norm = A / alpha

    W_q, S_vec, Vh_q = svd(A_norm)
    V_q = Vh_q.conj().T
    s_clamped = np.clip(S_vec, 0.0, 1.0)

    K_pad = 1 << int(np.ceil(np.log2(max(K, 1))))
    K_pad = max(1, min(K_pad, N))

    D_qiskit = np.ones(N, dtype=complex)
    for k in range(K):
        D_qiskit[k] = np.exp(1j * np.arccos(s_clamped[k]))
    for k in range(K, N):
        D_qiskit[k] = 1j

    iso_V = _isometry_first_k_columns(V_q, K_pad, n)

    # qiskit's ``Isometry`` synthesises the first ``K_pad`` columns
    # *exactly* (no residual per-column phase in qiskit >= 1.0), so the
    # diagonal residual ``Delta`` between ``iso_V`` and ``V_q[:, :K_pad]``
    # is identity to numerical precision.  Earlier revisions verified
    # this empirically via ``Operator(iso_V).data``, but that costs
    # ``O(4^n)`` time and memory; for n = 8 this single call dominated
    # the total run-time of ``siable_low_rank``.  We skip the check and
    # set ``Delta = I`` directly.  If a future qiskit release reintroduces
    # per-column phases this assertion will be caught by the block-
    # encoding correctness tests (Test 3, Test 5).
    iso_W = _isometry_first_k_columns(W_q, K_pad, n)

    qc = QuantumCircuit(n + 1)
    qc.h(0)
    qc.compose(iso_V.inverse(), qubits=list(range(1, n + 1)), inplace=True)

    P = _bit_reverse_indices(n)
    D_for_apply = D_qiskit[P]
    D_mat = np.diag(D_for_apply)
    _apply_ucrz(qc, D_mat, n_local=n + 1, mode="M", base=0)

    qc.compose(iso_W, qubits=list(range(1, n + 1)), inplace=True)
    qc.h(0)

    Sigma_K = np.zeros(N)
    Sigma_K[:K] = s_clamped[:K]
    A_K_qiskit = W_q @ np.diag(Sigma_K) @ Vh_q
    _fix_global_phase(qc, A_K_qiskit, n)

    return qc, alpha


# ---------------------------------------------------------------------------
# Gate wrappers (qclib-style)
# ---------------------------------------------------------------------------

class SiableBlockEncoding(Gate):
    """
    qclib-style block-encoding gate using SIABLE (full-rank).

    Example
    -------
    >>> from siable import SiableBlockEncoding
    >>> import numpy as np
    >>> A = np.random.randn(8, 8) + 1j * np.random.randn(8, 8)
    >>> gate = SiableBlockEncoding(A)
    >>> # `gate.alpha` is the subnormalisation factor.
    >>> # `gate.definition` is the corresponding QuantumCircuit.
    """

    def __init__(self, matrix, label: str = "SIABLE"):
        A = np.asarray(matrix, dtype=complex)
        n = int(round(np.log2(A.shape[0])))
        super().__init__("siable", n + 1, [], label=label)
        self._matrix = A
        qc, alpha = siable(A)
        self._qc = qc
        self.alpha = alpha

    def _define(self):
        self.definition = self._qc


class SiableLowRankBlockEncoding(Gate):
    """
    qclib-style block-encoding gate using SIABLE-low-rank (Theorem 4).

    Parameters
    ----------
    matrix : (N, N) array_like
    K : int
        Target rank.
    """

    def __init__(self, matrix, K: int, label: str = "SIABLE-LR"):
        A = np.asarray(matrix, dtype=complex)
        n = int(round(np.log2(A.shape[0])))
        super().__init__("siable_lr", n + 1, [], label=label)
        self._matrix = A
        self._K = K
        qc, alpha = siable_low_rank(A, K)
        self._qc = qc
        self.alpha = alpha

    def _define(self):
        self.definition = self._qc
