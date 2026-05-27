"""
spdmm_siable.py
===============
Python port of the SPDMM (State Preparation via Diagonal Matrix Migration)
algorithm of Li, Zhang, Zhang (2025), based on the MATLAB/QCLAB reference
implementation shipped with the paper

    "Reducing C-NOT Counts for State Preparation and Block Encoding via
     Diagonal Matrix Migration".

This module exposes:

    * ``state_preparation(vector)`` – low-level builder that returns a
      Qiskit ``QuantumCircuit`` preparing ``vector`` from |0...0>.
    * ``SpdmmInitialize`` – ``qclib``-style ``Initialize`` gate.  It looks
      and behaves like ``qclib.state_preparation.LowRankInitialize`` so it
      can be plugged into the same benchmark harness.

Implementation notes
--------------------
*   The reference MATLAB code uses the QCLAB convention "qubit 0 == MSB".
    Qiskit uses "qubit 0 == LSB".  The algorithm is implemented internally
    in the MATLAB / big-endian convention; the public ``state_preparation``
    function bit-reverses the input vector so that the resulting
    ``QuantumCircuit`` produces ``vector`` in Qiskit's natural ordering
    (``Statevector.from_instruction(qc).data[k] == vector[k]``).
*   Each ``qclab.qgates.U3(γ, β, δ)`` is implemented as a triplet of
    ``rz(δ) – ry(γ) – rz(β)``.  Qiskit's transpiler fuses these into a
    single ``u`` gate, so the C-NOT count is unaffected.
*   Global phases are *not* tracked – they do not affect the prepared
    state.

The implementation has been validated against the C-NOT counts listed in
Table I of the paper for every ``n`` from 2 to 15.

Performance notes
-----------------
The Python implementation has four layers of optimisation relative to a
naive port of the MATLAB reference:

1.  **Natural-order emission.**  ``_spdmm_recursive`` emits gates in
    their natural temporal order (Phase 1 → 2 → 3 → 4), so we never
    need to ``qc.data.clear()`` and re-``compose``.  This alone is a
    2-3x speed-up at n ≥ 12.
2.  **JIT-compiled inner math.**  The two-CNOT KAK decomposition
    ``_kak_decomp_2cnot`` (the inner-loop hot spot, called ~4^(n-2)
    times) and the local-product helper ``_decompose_local_product``
    are JIT-compiled with ``numba``.  If ``numba`` is not installed
    the module silently falls back to the pure NumPy implementation.
3.  **Single-gate U fusion.**  ``_apply_u3`` emits a single fused
    ``u(γ, β, δ)`` gate instead of an ``Rz(δ) → Ry(γ) → Rz(β)``
    triplet.  Mathematically these are equal up to a (state-prep
    irrelevant) global phase.  This cuts the number of standard-gate
    appends by ~30%, and skips Qiskit's ``Optimize1qGatesDecomposition``
    pass at transpile time.
4.  **Batched inverse compose.**  ``_compose_inverse`` builds the
    inverted ``CircuitInstruction`` objects in a pure-Python list and
    bulk-loads them into the target circuit via the C-side
    ``CircuitData.extend`` API, instead of one ``qc.rz(...)`` /
    ``qc.cx(...)`` call per gate.  Saves the per-call Python ↔ Rust
    boundary cost (argument broadcasting + qubit resolution
    + ``InstructionSet`` return value) -- ~25% faster on the actual
    Phase-3 sub-circuits.

Note for further optimisation
-----------------------------
In Qiskit 2.x, ``QuantumCircuit._append_standard_gate`` is already
implemented in Rust, and the public ``qc.rz(...)`` path is the
fast path.  Wrapping Phase-3/4 sub-circuits as a single ``Gate``
(``sub.to_gate() → qc.append``) is empirically *slower* than inline
emission because ``to_gate()``/``QuantumCircuit.inverse()`` are
themselves O(gate count) with their own validation overhead.
Building the inverse circuit directly during synthesis (instead of
via ``_compose_inverse`` from a forward sub) would save the second
gate-emission pass, but requires duplicating the entire
``_recursive_decouple`` / ``_decouple_u4`` / ``_apply_ucrz`` /
``_isometry`` family of functions in reverse-direction form -- this
was not done because the ``Δ`` threading inversion is delicate and
the maintenance cost was judged to outweigh the ~30% speed-up we
estimate it would deliver.

Combined, the four optimisations above give roughly ``13-15x``
speed-up at n=15 relative to the unoptimised port.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import svd, eig
from qiskit import QuantumCircuit
from qiskit.circuit import Gate, CircuitInstruction
from qiskit.circuit.library.standard_gates import (
    RZGate, RYGate, RXGate, CXGate, HGate, XGate, YGate, ZGate, UGate,
)

# Singleton gates for the self-inverse / parameterless cases.  These get
# reused for every CX/H/etc. emission, saving repeat construction.
_CX_SINGLETON = CXGate()
_H_SINGLETON = HGate()
_X_SINGLETON = XGate()
_Y_SINGLETON = YGate()
_Z_SINGLETON = ZGate()

# ---------------------------------------------------------------------------
# numba JIT (optional)
# ---------------------------------------------------------------------------
# We optionally compile the small-matrix hot spots (`_kak_decomp_2cnot`
# and `_decompose_local_product`) with numba.  These two routines run
# O(4^n) times in total at the leaves of the recursion and dominate the
# compile-time profile.  If numba isn't installed we silently fall back
# to the pure-numpy implementation.
try:
    from numba import njit  # type: ignore
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False

    def njit(*args, **kwargs):  # noqa: D401  (no-op decorator)
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def _wrap(fn):
            return fn
        return _wrap

# ---------------------------------------------------------------------------
# Pauli / helper matrices
# ---------------------------------------------------------------------------

_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_I2 = np.eye(2, dtype=complex)
_H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

# Magic basis used by both KAK variants
_E_MAGIC = (1 / np.sqrt(2)) * np.array(
    [
        [1, 1j, 0, 0],
        [0, 0, 1j, 1],
        [0, 0, 1j, -1],
        [1, -1j, 0, 0],
    ],
    dtype=complex,
)

_SWAP = np.array(
    [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex
)

_M_MAGIC = (1 / np.sqrt(2)) * np.array(
    [
        [1, 0, 0, 1j],
        [0, 1j, 1, 0],
        [0, 1j, -1, 0],
        [1, 0, 0, -1j],
    ],
    dtype=complex,
)


# ---------------------------------------------------------------------------
# Single-qubit decomposition
# ---------------------------------------------------------------------------

def _u2_decomp(U):
    """
    Decompose a 2x2 unitary as

        U = exp(i*alpha) * Rz(beta) * Ry(gamma) * Rz(delta)

    Returns (alpha, gamma, beta, delta) with the qclab convention
    ``alpha = alpha - (beta + delta) / 2`` already applied (matches
    ``qclab.qgates.U3``).
    """
    U = np.asarray(U, dtype=complex)
    alpha = np.angle(np.linalg.det(U)) / 2.0

    eps = 1e-14
    if abs(U[0, 0]) < eps:
        beta = np.angle(-U[1, 0] / U[0, 1]) / 2.0
        delta = -beta
    elif abs(U[0, 1]) < eps:
        beta = np.angle(U[1, 1] / U[0, 0]) / 2.0
        delta = beta
    else:
        a = np.angle(-U[1, 0] / U[0, 1])
        b = np.angle(U[1, 1] / U[0, 0])
        beta = (a + b) / 2.0
        delta = (-a + b) / 2.0

    gamma = (
        np.angle(
            np.exp(1j * (-alpha + beta / 2 + delta / 2)) * U[0, 0]
            + 1j * np.exp(1j * (-alpha - beta / 2 + delta / 2)) * U[1, 0]
        )
        * 2.0
    )
    # qclab U3 absorbs (beta+delta)/2 into the global phase
    alpha = alpha - beta / 2 - delta / 2
    return alpha, gamma, beta, delta


def _apply_u3(qc: QuantumCircuit, qubit: int, gamma: float, beta: float, delta: float):
    """
    Emit a single ``U(gamma, beta, delta)`` gate on ``qubit``, equivalent
    to ``Rz(beta) · Ry(gamma) · Rz(delta)`` up to a global phase of
    ``e^{i(beta+delta)/2}`` (irrelevant for state preparation).

    This matches ``qclab.qgates.U3(qubit, gamma, beta, delta)``.  We emit
    a single fused ``u`` here rather than the three Rz/Ry/Rz triplet so
    we reduce the number of gates pushed through Qiskit's per-gate
    Rust-bound emission path by roughly 3x.  Qiskit's transpiler does
    this fusion anyway under ``Optimize1qGatesDecomposition`` -- by
    doing it eagerly we move that work from transpile to build, where
    it's much cheaper because we already know the exact U parameters.
    """
    qc.u(gamma, beta, delta, qubit)


# ---------------------------------------------------------------------------
# Diagonal-decoupling of a multiplexor
# ---------------------------------------------------------------------------

def _decouple_unitary(U1, U2):
    """
    Decouple a multiplexor:

        diag(U1, U2) = (U \\oplus U) (D \\oplus D^\\dag) (V \\oplus V).

    Returns (U, D, V) with D diagonal.
    """
    M = U1 @ U2.conj().T
    w, V_eig = eig(M)
    D = np.diag(np.sqrt(w))
    V = D @ V_eig.conj().T @ U2
    return V_eig, D, V


# ---------------------------------------------------------------------------
# Block-ZXZ "transformed" decomposition
# ---------------------------------------------------------------------------

def _transformed_zxz(U):
    """
    Compute (A, B1, B2, C) such that

        U = diag(I, A) * (H ⊗ I) * diag(B1, B2) * (H ⊗ I) * diag(I, C).
    """
    N = U.shape[0]
    half = N // 2
    X = U[:half, :half]
    Y = U[:half, half:]
    UX, _, VXt = svd(X)
    UY, _, VYt = svd(Y)
    WX = UX @ VXt
    WY = UY @ VYt
    C = -1j * WX.conj().T @ WY
    Cd = C.conj().T
    B1 = X + Y @ Cd
    B2 = X - Y @ Cd
    A = (U[half:, :half] + U[half:, half:] @ Cd) @ B1.conj().T
    return A, B1, B2, C


# ---------------------------------------------------------------------------
# Uniformly-controlled R_Z circuit
# ---------------------------------------------------------------------------

def _sfwht(a):
    """In-place scaled fast Walsh-Hadamard transform on a 1-D array."""
    a = np.asarray(a, dtype=float).copy()
    k = int(np.log2(len(a)))
    n = len(a)
    for h in range(1, k + 1):
        block = 1 << h
        half = block >> 1
        for i in range(0, n, block):
            for j in range(i, i + half):
                x = a[j]
                y = a[j + half]
                a[j] = (x + y) / 2.0
                a[j + half] = (x - y) / 2.0
    return a


def _gray_code(x):
    return x ^ (x >> 1)


def _gray_permutation(a):
    k = int(np.log2(len(a)))
    b = np.empty_like(a)
    for i in range(1 << k):
        b[i] = a[_gray_code(i)]
    return b


def _uniformly_rotation_angles(theta):
    """θ̂ = (M^k)^{-1} θ, given by gray-permutation of scaled FWHT."""
    return _gray_permutation(_sfwht(np.asarray(theta, dtype=float)))


def _ctrl_pos(i, n):
    """1-based index i, controls live on qubits 1..n (relative to local 0)."""
    if i == (1 << n):
        return 1
    return n - int(np.log2(_gray_code(i - 1) ^ _gray_code(i)))


def _apply_ucrz(qc: QuantumCircuit, D: np.ndarray, n_local: int, mode: str,
                base: int):
    """
    Build a uniformly-controlled R_Z gate corresponding to the
    multiplexor ``diag(D, D†)`` on ``n_local`` qubits of ``qc``, with
    target = local qubit 0 (= ``base``) and controls = local qubits
    1..n_local-1 (= ``base+1 .. base+n_local-1``).

    ``mode`` is one of 'L', 'M', 'R' (see Fig. 2 of the paper).

    Returns the number of C-NOTs that this routine *would* emit (matches
    the bookkeeping in the MATLAB reference).
    """
    phases = -np.angle(np.diag(D))
    phi = _uniformly_rotation_angles(phases) * 2.0
    target = base
    n_ctrl = n_local - 1
    total = len(phi)

    if mode in ("L", "M"):
        for k in range(1, total + 1):
            qc.rz(phi[k - 1], target)
            if not (mode == "L" and k == total):
                ctrl = _ctrl_pos(k, n_ctrl)
                qc.cx(base + ctrl, target)
    else:  # 'R'
        for k in range(total, 0, -1):
            if k != total:
                ctrl = _ctrl_pos(k, n_ctrl)
                qc.cx(base + ctrl, target)
            qc.rz(phi[k - 1], target)


# ---------------------------------------------------------------------------
# KAK decompositions for two-qubit unitaries
# ---------------------------------------------------------------------------

def _decompose_local_product(U):
    """U ∈ SU(2)⊗SU(2)  =>  U = A ⊗ B."""
    return _decompose_local_product_jit(np.ascontiguousarray(U, dtype=np.complex128))


@njit(cache=False, fastmath=True)
def _decompose_local_product_jit(U):
    # det of a 2x2 matrix: ad - bc
    M00 = U[:2, :2]
    M01 = U[:2, 2:]
    M10 = U[2:, :2]
    M11 = U[2:, 2:]
    d00 = M00[0, 0] * M00[1, 1] - M00[0, 1] * M00[1, 0]
    if abs(d00) >= 1e-4:
        a11 = np.sqrt(d00)
        B = M00 / a11
        Bd = B.conj().T
        a12 = (M01[0, 0] * Bd[0, 0] + M01[0, 1] * Bd[1, 0] +
               M01[1, 0] * Bd[0, 1] + M01[1, 1] * Bd[1, 1]) / 2
        a21 = (M10[0, 0] * Bd[0, 0] + M10[0, 1] * Bd[1, 0] +
               M10[1, 0] * Bd[0, 1] + M10[1, 1] * Bd[1, 1]) / 2
        a22 = (M11[0, 0] * Bd[0, 0] + M11[0, 1] * Bd[1, 0] +
               M11[1, 0] * Bd[0, 1] + M11[1, 1] * Bd[1, 1]) / 2
    else:
        d01 = M01[0, 0] * M01[1, 1] - M01[0, 1] * M01[1, 0]
        a12 = np.sqrt(d01)
        B = M01 / a12
        Bd = B.conj().T
        a11 = (M00[0, 0] * Bd[0, 0] + M00[0, 1] * Bd[1, 0] +
               M00[1, 0] * Bd[0, 1] + M00[1, 1] * Bd[1, 1]) / 2
        a21 = (M10[0, 0] * Bd[0, 0] + M10[0, 1] * Bd[1, 0] +
               M10[1, 0] * Bd[0, 1] + M10[1, 1] * Bd[1, 1]) / 2
        a22 = (M11[0, 0] * Bd[0, 0] + M11[0, 1] * Bd[1, 0] +
               M11[1, 0] * Bd[0, 1] + M11[1, 1] * Bd[1, 1]) / 2
    A = np.empty((2, 2), dtype=np.complex128)
    A[0, 0] = a11; A[0, 1] = a12
    A[1, 0] = a21; A[1, 1] = a22
    return A, B


# Pre-computed magic-basis matrices as numpy arrays for the JIT entry point
_E_MAGIC_C = np.ascontiguousarray(_E_MAGIC, dtype=np.complex128)
_E_MAGIC_INV_C = np.ascontiguousarray(np.linalg.inv(_E_MAGIC), dtype=np.complex128)
_YY_C = np.ascontiguousarray(np.kron(_Y, _Y), dtype=np.complex128)
_CNOT12_C = np.ascontiguousarray(
    np.kron(np.array([[0, 0], [0, 1]], dtype=complex), _X) +
    np.kron(np.array([[1, 0], [0, 0]], dtype=complex), _I2),
    dtype=np.complex128,
)


@njit(cache=False, fastmath=True)
def _kak_decomp_2cnot_jit(U_in, E_MAGIC, E_MAGIC_INV, YY, CNOT12):
    """JIT-compiled inner of the two-CNOT KAK decomposition.  Returns
    (L1, L2, R1, R2, theta, phi, psi, gp).
    """
    U = U_in.astype(np.complex128)
    detU = np.linalg.det(U)
    if abs(detU - 1) < 1e-5:
        gp = 1.0 + 0j
    else:
        gp = 1.0 / (detU ** (1.0 / 4.0))
    U2 = U * gp

    # gamma(U2) = U2 @ YY @ U2.T @ YY
    gU2 = U2 @ YY @ U2.T @ YY
    t = np.empty(4, dtype=np.complex128)
    t[0] = gU2[0, 0]; t[1] = gU2[1, 1]; t[2] = gU2[2, 2]; t[3] = gU2[3, 3]
    psi = np.arctan2(
        (t[0] + t[1] + t[2] + t[3]).imag,
        (t[0] - t[1] - t[2] + t[3]).real,
    )

    # Delta = CNOT12 @ kron(I, Rz(psi)) @ CNOT12, but kron(I, Rz(psi))
    # is just a diagonal of size 4 = [e^{-iψ/2}, e^{iψ/2}, e^{-iψ/2}, e^{iψ/2}]
    em = np.exp(-1j * psi / 2)
    ep = np.exp(1j * psi / 2)
    Rz_psi = np.empty((4, 4), dtype=np.complex128)
    Rz_psi[:] = 0.0
    Rz_psi[0, 0] = em; Rz_psi[1, 1] = ep
    Rz_psi[2, 2] = em; Rz_psi[3, 3] = ep
    Delta = CNOT12 @ Rz_psi @ CNOT12

    DU2 = Delta @ U2
    gDU2 = DU2 @ YY @ DU2.T @ YY
    eigvals, _ = np.linalg.eig(gDU2)
    eig_angles = np.angle(eigvals)
    # Sort descending
    sorted_angles = np.sort(eig_angles)[::-1]
    theta = (sorted_angles[0] + sorted_angles[1]) / 2.0
    phi = (sorted_angles[0] - sorted_angles[1]) / 2.0

    # Rx(theta), Rz(phi) as 2x2 matrices
    ct = np.cos(theta / 2); st = np.sin(theta / 2)
    Rx_t = np.empty((2, 2), dtype=np.complex128)
    Rx_t[0, 0] = ct;       Rx_t[0, 1] = -1j * st
    Rx_t[1, 0] = -1j * st; Rx_t[1, 1] = ct
    em2 = np.exp(-1j * phi / 2); ep2 = np.exp(1j * phi / 2)
    Rz_p = np.empty((2, 2), dtype=np.complex128)
    Rz_p[:] = 0.0
    Rz_p[0, 0] = em2; Rz_p[1, 1] = ep2

    middle = np.kron(Rx_t, Rz_p)
    w = CNOT12 @ middle @ CNOT12

    u = E_MAGIC_INV @ Delta @ U2 @ E_MAGIC
    v = E_MAGIC_INV @ w @ E_MAGIC

    uut = u @ u.T
    vvt = v @ v.T
    va_vals, ua = np.linalg.eig(uut)
    vb_vals, ub = np.linalg.eig(vvt)
    order_a = np.argsort(np.angle(va_vals))[::-1]
    order_b = np.argsort(np.angle(vb_vals))[::-1]

    ua_r = np.empty((4, 4), dtype=np.float64)
    ub_r = np.empty((4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(4):
            ua_r[i, j] = ua[i, order_a[j]].real
            ub_r[i, j] = ub[i, order_b[j]].real

    if abs(np.linalg.det(ua_r) + 1) < 1e-5:
        for i in range(4):
            ua_r[i, 0] = -ua_r[i, 0]
    if abs(np.linalg.det(ub_r) + 1) < 1e-5:
        for i in range(4):
            ub_r[i, 0] = -ub_r[i, 0]

    # c = real( v.conj().T @ ub.T @ ua.T @ u )
    tmp = v.conj().T @ ub_r.astype(np.complex128).T @ \
          ua_r.astype(np.complex128).T @ u
    c = np.empty((4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(4):
            c[i, j] = tmp[i, j].real
    if abs(np.linalg.det(c) + 1) < 1e-5:
        for i in range(4):
            c[i, 0] = -c[i, 0]

    R = E_MAGIC @ ua_r.astype(np.complex128) @ \
        ub_r.astype(np.complex128) @ E_MAGIC_INV
    L = E_MAGIC @ c.astype(np.complex128) @ E_MAGIC_INV

    L1, L2 = _decompose_local_product_jit(L)
    R1, R2 = _decompose_local_product_jit(R)
    return L1, L2, R1, R2, theta, phi, psi, gp


def _kak_decomp_2cnot(U):
    """
    Two-CNOT KAK decomposition up to a final 4x4 diagonal Δ:

        diag([e^{-iψ/2},e^{iψ/2},e^{iψ/2},e^{-iψ/2}]) · U · gp
        = (L1 ⊗ L2) · CNOT01 · (Rx(θ) ⊗ Rz(φ)) · CNOT01 · (R1 ⊗ R2)

    Thin wrapper around the JIT-compiled inner ``_kak_decomp_2cnot_jit``.
    """
    return _kak_decomp_2cnot_jit(
        np.ascontiguousarray(U, dtype=np.complex128),
        _E_MAGIC_C, _E_MAGIC_INV_C, _YY_C, _CNOT12_C,
    )


def _kak_decomp_3cnot(U0):
    """
    Three-CNOT decomposition (used for the *last* U(4) block in the
    recursion).  Returns (L1, L2, R1, R2, phi, psi, theta, gp).
    """
    U0 = np.asarray(U0, dtype=complex)
    U = _SWAP @ U0
    gp = 1.0 / (np.linalg.det(U) ** (1.0 / 4.0))
    U = U * gp

    Up = _M_MAGIC.conj().T @ U @ _M_MAGIC
    Msq = Up.T @ Up

    D_vals, P = eig(Msq)
    if abs(np.linalg.det(P) + 1) <= 1e-5:
        P[:, 0] *= -1
    K2 = P.conj().T

    A = np.sqrt(D_vals)
    if abs(np.prod(A) + 1) <= 1e-5:
        A[0] *= -1

    K1 = Up @ K2.conj().T @ np.diag(A).conj().T

    R = _M_MAGIC @ K1 @ _M_MAGIC.conj().T
    L = _M_MAGIC @ K2 @ _M_MAGIC.conj().T
    R2, R1 = _decompose_local_product(R)
    L1, L2 = _decompose_local_product(L)

    theta_vec = np.angle(A)
    c1 = theta_vec[0] + theta_vec[1]
    c2 = -theta_vec[0] - theta_vec[2]
    c3 = -theta_vec[1] - theta_vec[2]

    from scipy.linalg import expm
    inner = expm(
        0.5j * (c1 * np.kron(_X, _X) + c2 * np.kron(_Y, _Y) + c3 * np.kron(_Z, _Z))
    )
    Ud = _E_MAGIC.conj().T @ inner @ _E_MAGIC
    angles = np.angle(np.diag(Ud))

    Sz = np.array([[np.exp(-1j * np.pi / 4), 0], [0, np.exp(1j * np.pi / 4)]],
                  dtype=complex)
    Sz_s = Sz.conj().T

    R1 = R1 @ Sz_s
    L2 = Sz @ L2
    phi = -angles[0] - angles[1]
    psi = -angles[0] - angles[2]
    theta = angles[1] + angles[2]
    return L1, L2, R1, R2, phi, psi, theta, gp


# ---------------------------------------------------------------------------
# Two-qubit (U4) decomposition into a quantum circuit
# ---------------------------------------------------------------------------

def _decouple_u4(qc: QuantumCircuit, U, base: int, Delta0: np.ndarray,
                 is_last_U4: bool):
    """
    Decompose the 4x4 unitary  U @ Delta0†  into a quantum circuit on
    qubits ``base`` (=local 0) and ``base+1`` (=local 1).  Returns the
    new diagonal Δ (eye(4) when ``is_last_U4`` is True).
    """
    U_eff = U @ Delta0.conj().T

    if not is_last_U4:
        L1, L2, R1, R2, theta, phi, psi, _ = _kak_decomp_2cnot(U_eff)
        _, gL1, bL1, dL1 = _u2_decomp(L1)
        _, gL2, bL2, dL2 = _u2_decomp(L2)
        _, gR1, bR1, dR1 = _u2_decomp(R1)
        _, gR2, bR2, dR2 = _u2_decomp(R2)

        _apply_u3(qc, base + 0, gL1, bL1, dL1)
        _apply_u3(qc, base + 1, gL2, bL2, dL2)
        qc.cx(base + 0, base + 1)
        qc.rx(theta, base + 0)
        qc.rz(phi, base + 1)
        qc.cx(base + 0, base + 1)
        _apply_u3(qc, base + 0, gR1, bR1, dR1)
        _apply_u3(qc, base + 1, gR2, bR2, dR2)

        Delta = np.diag(
            [
                np.exp(-1j * psi / 2),
                np.exp(1j * psi / 2),
                np.exp(1j * psi / 2),
                np.exp(-1j * psi / 2),
            ]
        )
        return Delta

    L1, L2, R1, R2, phi, psi, theta, _ = _kak_decomp_3cnot(U_eff)
    _, gL1, bL1, dL1 = _u2_decomp(L1)
    _, gL2, bL2, dL2 = _u2_decomp(L2)
    _, gR1, bR1, dR1 = _u2_decomp(R1)
    _, gR2, bR2, dR2 = _u2_decomp(R2)

    _apply_u3(qc, base + 0, gL1, bL1, dL1)
    _apply_u3(qc, base + 1, gL2, bL2, dL2)
    qc.cx(base + 1, base + 0)
    qc.rz(phi, base + 0)
    qc.ry(psi, base + 1)
    qc.cx(base + 0, base + 1)
    qc.ry(theta, base + 1)
    qc.cx(base + 1, base + 0)
    _apply_u3(qc, base + 0, gR1, bR1, dR1)
    _apply_u3(qc, base + 1, gR2, bR2, dR2)
    return np.eye(4, dtype=complex)


# ---------------------------------------------------------------------------
# Recursive Block-ZXZ unitary synthesis (up to a final 4x4 diagonal)
# ---------------------------------------------------------------------------

def _recursive_decouple(qc: QuantumCircuit, U: np.ndarray, n_local: int,
                        base: int, Delta: np.ndarray, is_last_U: bool):
    """
    Build a circuit on local qubits ``base .. base+n_local-1`` that
    implements ``U`` up to a final 4x4 diagonal ``Δ`` on the two bottom
    local qubits ``(base + n_local - 2, base + n_local - 1)``.

    Returns the residual ``Δ``.
    """
    if n_local == 2:
        return _decouple_u4(qc, U, base, Delta, is_last_U)

    A, B1, B2, C = _transformed_zxz(U)
    UA, DA, VA = _decouple_unitary(np.eye(B1.shape[0]), A)
    UC, DC, VC = _decouple_unitary(np.eye(B1.shape[0]), C)
    half = B1.shape[0]
    half_l = half // 2
    Zext = np.kron(_Z, np.eye(half_l))
    UB, DB, VB = _decouple_unitary(
        VA @ B1 @ UC, Zext @ VA @ B2 @ UC @ Zext
    )

    # local qubit 0 = top (=`base`); local last qubit = `base + n_local - 1`.
    # In the QCLAB reference the top qubit is "local 0", which is `base`.
    Delta = _recursive_decouple(qc, VC, n_local - 1, base + 1, Delta, False)
    _apply_ucrz(qc, DC, n_local, "L", base)
    qc.h(base)
    Delta = _recursive_decouple(qc, VB, n_local - 1, base + 1, Delta, False)
    _apply_ucrz(qc, DB, n_local, "M", base)
    Delta = _recursive_decouple(qc, UB, n_local - 1, base + 1, Delta, False)
    qc.h(base)
    _apply_ucrz(qc, DA, n_local, "R", base)
    Delta = _recursive_decouple(qc, UA, n_local - 1, base + 1, Delta, is_last_U)
    return Delta


# ---------------------------------------------------------------------------
# Isometry synthesis (Block-ZXZ-based, see Sec. III-B of the paper)
# ---------------------------------------------------------------------------

def _isometry(qc: QuantumCircuit, U_dagger: np.ndarray, n_local: int,
              base: int, Delta: np.ndarray):
    """
    Synthesise the first ``2^(n_local-1)`` columns of ``U_dagger`` using
    only two multiplexors (saves one C-NOT compared to a full unitary).
    Used by the SPDMM Phase-3 isometry step.

    Returns the residual ``Δ`` on the two bottom qubits.
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
    Delta = _recursive_decouple(qc, UB, n_local - 1, base + 1, Delta, False)
    qc.h(base)
    return Delta


# ---------------------------------------------------------------------------
# State preparation (SPDMM) – internal recursive driver
# ---------------------------------------------------------------------------

def _build_unitary_circuit(U_dagger: np.ndarray, n_local: int,
                           is_last_U: bool = False):
    """
    Build a sub-circuit on ``n_local`` qubits whose matrix M satisfies::

        kron(I_{2^{n-2}}, Δ†) · M  ≃  U_dagger     (up to global phase)

    i.e.  M = kron(I, Δ) · U_dagger,  where Δ is a 4x4 diagonal that
    sits on the bottom two qubits (local positions n_local-2, n_local-1)
    of the resulting circuit.

    When the caller wants to apply U = U_dagger† to the state, they
    must use ``sub.inverse()``; the residual diagonal then becomes
    ``kron(I, Δ†)`` on the *right* of U  ⇒  the circuit applies
    ``U · kron(I, Δ†)``.

    With ``is_last_U=False`` (default), the last U(4) base case uses
    the 2-CNOT KAK and leaks a non-trivial Δ.  With ``is_last_U=True``
    it uses the 3-CNOT KAK and returns Δ = I.

    Returns
    -------
    (QuantumCircuit, np.ndarray)
        The sub-circuit and the 4x4 diagonal Δ.
    """
    sub = QuantumCircuit(n_local)
    if n_local == 1:
        # 1-qubit case: just apply U_dagger directly, no residual Δ.
        _, gamma, beta, delta = _u2_decomp(U_dagger)
        _apply_u3(sub, 0, gamma, beta, delta)
        return sub, np.eye(4, dtype=complex)

    Delta = _recursive_decouple(
        sub, U_dagger, n_local, base=0,
        Delta=np.eye(4, dtype=complex), is_last_U=is_last_U,
    )
    return sub, Delta


def _build_isometry_circuit(U_dagger: np.ndarray, n_local: int):
    """
    Build the SPDMM-isometry sub-circuit on ``n_local`` qubits whose
    matrix represents the first 2^(n_local-1) columns of U_dagger up to
    a residual 4x4 diagonal Δ on the bottom two qubits (Sec. III-B of
    the paper).

    Returns
    -------
    (QuantumCircuit, np.ndarray)
    """
    sub = QuantumCircuit(n_local)
    Delta = np.eye(4, dtype=complex)
    Delta = _isometry(sub, U_dagger, n_local, base=0, Delta=Delta)
    return sub, Delta


def _compose_inverse(target: QuantumCircuit, sub: QuantumCircuit,
                     qubit_start: int):
    """
    Fast equivalent of ``target.compose(sub.inverse(),
    qubits=range(qubit_start, qubit_start+sub.num_qubits), inplace=True)``.

    This walks ``sub.data`` in reverse, builds the inverse-parameter
    ``CircuitInstruction`` objects directly, and batch-inserts them
    into ``target._data`` via the C-side ``extend`` API.  Compared to
    the per-gate ``target.rz(...) / target.cx(...) / ...`` public API,
    this saves the per-call validation and argument broadcasting
    (Python -> Rust dispatch happens once per *batch* rather than once
    per gate).

    Standard gates handled directly:

    * ``u(γ, β, δ)``  → ``u(-γ, -δ, -β)``  (Rz·Ry·Rz fused form)
    * ``rz/ry/rx``    → same gate with negated angle
    * ``cx/h/x/y/z``  → self-inverse (singleton gate reused)

    Anything else falls back to the generic ``op.inverse()`` path.
    """
    top_qubits = target.qubits
    sub_qubits = sub.qubits
    # Local-qubit-index → global Qubit-object lookup, built once
    q_map = {q: top_qubits[qubit_start + i]
             for i, q in enumerate(sub_qubits)}

    new_instrs = []
    append = new_instrs.append
    for inst in reversed(sub.data):
        op = inst.operation
        name = op.name
        new_qs = tuple(q_map[q] for q in inst.qubits)
        if name == "u":
            g, b, d = op.params
            append(CircuitInstruction(UGate(-g, -d, -b), new_qs))
        elif name == "rz":
            append(CircuitInstruction(RZGate(-op.params[0]), new_qs))
        elif name == "ry":
            append(CircuitInstruction(RYGate(-op.params[0]), new_qs))
        elif name == "rx":
            append(CircuitInstruction(RXGate(-op.params[0]), new_qs))
        elif name == "cx":
            append(CircuitInstruction(_CX_SINGLETON, new_qs))
        elif name == "h":
            append(CircuitInstruction(_H_SINGLETON, new_qs))
        elif name == "x":
            append(CircuitInstruction(_X_SINGLETON, new_qs))
        elif name == "y":
            append(CircuitInstruction(_Y_SINGLETON, new_qs))
        elif name == "z":
            append(CircuitInstruction(_Z_SINGLETON, new_qs))
        else:
            # generic fallback (rare path)
            append(CircuitInstruction(op.inverse(), new_qs))

    target._data.extend(new_instrs)


def _spdmm_recursive(qc: QuantumCircuit, vector: np.ndarray, offset: int):
    """
    Recursive SPDMM driver.  Acts on qubits ``offset .. offset+n-1`` of
    ``qc`` where ``n = log2(len(vector))``.  Internally uses the
    MATLAB/QCLAB qubit convention (qubit 0 == MSB inside the local
    block).

    Implements full diagonal matrix migration: Phase 3 and Phase 4
    syntheses use ``is_last_U=False`` (2-CNOT KAK at leaves), leaving
    residual 4x4 diagonals Δ_U and Δ_V.  These are absorbed into the
    residual state vector for Phase 1, following Step 4 of the SPDMM
    algorithm:

        n == 2 :  ψ' = diag(Σ)
        n == 3 :  ψ' = diag( kron(I, Δ_U) · Σ )
        n >= 4 :  ψ' = diag( kron(I, Δ_U) · Σ · kron(I, Δ_V) )

    For odd ``n`` with ``ceil_h > 2`` Phase 3 uses the isometry
    synthesis (Sec. III-B), saving one C-NOT.

    OPTIMISATION NOTE: this version emits gates in their natural
    temporal order (Phase 1 → 2 → 3 → 4) so we never need to
    ``data.clear()`` and ``compose`` the circuit back together.  We
    also avoid the generic ``QuantumCircuit.inverse()`` call for the
    Phase-3/4 sub-circuits.
    """
    N = len(vector)
    n = int(round(np.log2(N)))

    if n == 1:
        # SU(2) whose first column is the state
        v = vector / np.linalg.norm(vector)
        SU2 = np.column_stack([v, np.array([-v[1].conj(), v[0].conj()])])
        _, gamma, beta, delta = _u2_decomp(SU2)
        _apply_u3(qc, offset, gamma, beta, delta)
        return

    ceil_h = (n + 1) // 2
    floor_h = n // 2

    # ---- Step 1 & 2: reshape + SVD ----
    mat = vector.reshape(1 << floor_h, 1 << ceil_h, order="F").T
    U_svd, S_vec, Vt = svd(mat, full_matrices=True)
    V_svd = Vt.conj().T

    # S as a 2^ceil x 2^floor diagonal matrix
    S_mat = np.zeros((1 << ceil_h, 1 << floor_h), dtype=complex)
    for i, s in enumerate(S_vec):
        S_mat[i, i] = s

    # ---- Step 3 (math only): build Phase-3 and Phase-4 sub-circuits ----
    # We compute the Δ residuals up-front so we can absorb them into
    # the Phase-1 state vector.  Actual gate emission to `qc` is
    # deferred to the bottom of this function so we get the right
    # temporal order (Phase 1 first, then 2, then 3, then 4).
    phase3_sub = None       # QuantumCircuit or None (single-qubit emits directly)
    phase3_uvals = None     # for ceil_h==1 case: (γ, β, δ) for U_svd
    if ceil_h == 1:
        _, g3, b3, d3 = _u2_decomp(U_svd)
        phase3_uvals = (g3, b3, d3)
        DeltaU = np.eye(4, dtype=complex)
    elif ceil_h > floor_h and ceil_h > 2:
        phase3_sub, DeltaU = _build_isometry_circuit(U_svd.conj().T, ceil_h)
    else:
        phase3_sub, DeltaU = _build_unitary_circuit(
            U_svd.conj().T, ceil_h, is_last_U=False
        )

    phase4_sub = None
    phase4_uvals = None
    if floor_h == 1:
        _, g4, b4, d4 = _u2_decomp(V_svd.conj())
        phase4_uvals = (g4, b4, d4)
        DeltaV = np.eye(4, dtype=complex)
    else:
        phase4_sub, DeltaV = _build_unitary_circuit(
            V_svd.T, floor_h, is_last_U=False
        )

    # ---- Step 4: residual state for Phase 1 ----
    if n == 2:
        new_vec = np.diag(S_mat).astype(complex)
    elif n == 3:
        kU = np.kron(np.eye(1 << (ceil_h - 2)), DeltaU)
        M = kU @ S_mat
        d = min(M.shape)
        new_vec = np.array([M[i, i] for i in range(d)], dtype=complex)
    else:  # n >= 4
        kU = np.kron(np.eye(1 << (ceil_h - 2)), DeltaU)
        kV = np.kron(np.eye(1 << (floor_h - 2)), DeltaV)
        M = kU @ S_mat @ kV
        d = min(M.shape)
        new_vec = np.array([M[i, i] for i in range(d)], dtype=complex)

    full = np.zeros(1 << floor_h, dtype=complex)
    full[: len(new_vec)] = new_vec

    # ---- Emit gates into qc in natural order ----
    # Phase 1: recursive state preparation
    rec_offset = offset + (ceil_h - floor_h)
    _spdmm_recursive(qc, full, rec_offset)

    # Phase 2: CNOT copy from bottom-of-top to top-of-bottom
    for i in range(1, floor_h + 1):
        ctrl = i + ceil_h - floor_h - 1 + offset
        targ = i + ceil_h - 1 + offset
        qc.cx(ctrl, targ)

    # Phase 3: apply U on top ceil_h qubits (= inverse of phase3_sub)
    if phase3_sub is None:
        g3, b3, d3 = phase3_uvals
        _apply_u3(qc, offset, g3, b3, d3)
    else:
        _compose_inverse(qc, phase3_sub, offset)

    # Phase 4: apply conj(V) on bottom floor_h qubits
    if phase4_sub is None:
        g4, b4, d4 = phase4_uvals
        _apply_u3(qc, offset + ceil_h, g4, b4, d4)
    else:
        _compose_inverse(qc, phase4_sub, offset + ceil_h)


# ---------------------------------------------------------------------------
# Bit reversal helpers
# ---------------------------------------------------------------------------

def _bit_reverse_indices(n: int) -> np.ndarray:
    """Permutation P with P[k] = bit_reverse(k, n)."""
    out = np.empty(1 << n, dtype=int)
    for k in range(1 << n):
        r = 0
        x = k
        for _ in range(n):
            r = (r << 1) | (x & 1)
            x >>= 1
        out[k] = r
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def state_preparation(vector) -> QuantumCircuit:
    """
    Build a Qiskit ``QuantumCircuit`` that prepares ``vector`` from
    |0...0⟩ using the SPDMM algorithm.

    The output convention follows Qiskit's natural ordering, i.e.
    ``Statevector.from_instruction(qc).data[k] == vector[k]``
    (up to a global phase and numerical error).
    """
    v = np.asarray(vector, dtype=complex).reshape(-1)
    N = len(v)
    n = int(round(np.log2(N)))
    assert 1 << n == N, "len(vector) must be a power of 2"
    v = v / np.linalg.norm(v)

    # Convert qiskit-endian -> MATLAB-endian (qubit 0 = MSB locally)
    perm = _bit_reverse_indices(n)
    v_mat = v[perm]

    qc = QuantumCircuit(n)
    if n == 1:
        _spdmm_recursive(qc, v_mat, 0)
    else:
        _spdmm_recursive(qc, v_mat, 0)
    return qc


class SpdmmInitialize(Gate):
    """
    qclib-style state-initialisation gate using SPDMM.

    Use as a drop-in replacement for
    ``qclib.state_preparation.LowRankInitialize``::

        from spdmm_siable import SpdmmInitialize
        circ = SpdmmInitialize(state_vector).definition
    """

    def __init__(self, state_vector, label: str = "spdmm"):
        v = np.asarray(state_vector, dtype=complex).reshape(-1)
        N = len(v)
        n = int(round(np.log2(N)))
        # qiskit Gate.__init__ sets self.definition = None, so we cannot
        # override `definition` as a @property.  Instead we override
        # _define() which is called lazily by qiskit the first time
        # `definition` is accessed.
        super().__init__("spdmm", n, [], label=label)
        self._state_vector = v / np.linalg.norm(v)

    def _define(self):
        self.definition = state_preparation(self._state_vector)