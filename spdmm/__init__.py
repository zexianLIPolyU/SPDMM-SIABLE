"""
spdmm
=====
Reference Python implementation of two algorithms from
Li, Zhang, & Zhang (2025), *"Reducing C-NOT Counts for State Preparation
and Block Encoding via Diagonal Matrix Migration"* (IEEE TCAD).

    SPDMM   -- State Preparation via Diagonal Matrix Migration
    SIABLE  -- Single-Ancilla Block Encoding

Public API
----------
state_preparation(state)
    Prepare an n-qubit state via SPDMM. Returns a ``QuantumCircuit`` whose
    output state matches ``state`` up to a global phase.  Leading C-NOT
    constant: ``(11/12) * 2^n``.

siable(matrix)
    Block-encode an arbitrary ``2^n x 2^n`` matrix with a single ancilla
    qubit (qiskit qubit 0).  Returns ``(circuit, alpha)`` where the
    block-encoding relation is ``A == alpha * U[0::2, 0::2]``.
    Leading C-NOT constant: ``(11/48) * 4^(n+1)``.

siable_low_rank(matrix, K)
    Block-encode a rank-K matrix with a single ancilla, automatically
    selecting the cheapest of four sub-constructions (state-preparation,
    column-by-column isometry, Block-ZXZ half-isometry, or full-rank).

SiableBlockEncoding
    A ``qiskit.circuit.Gate`` wrapper exposing ``siable`` / ``siable_low_rank``
    as a composable circuit element.

SpdmmInitialize
    A ``qclib``-style adapter around ``state_preparation`` exposing the
    ``.definition`` attribute used by qclib's benchmark conventions.  Useful
    for plugging SPDMM into the same harness as ``qclib.state_preparation``
    methods like ``LowRankInitialize``.
"""

from qiskit import QuantumCircuit
import numpy as np

from .spdmm import state_preparation
from .siable import siable, siable_low_rank, SiableBlockEncoding

# The Rust hybrid (qsd-based leaf synthesis, +7 CNOTs at n=15, ~2x faster).
# Imported eagerly so ``from spdmm import SpdmmInitializeRust`` works; the
# benchmark runner also auto-detects the module by name.
from .spdmm_hybrid_rust import SpdmmInitializeRust, state_preparation_rust


# ---------------------------------------------------------------------------
# qclib-style adapter
# ---------------------------------------------------------------------------

class SpdmmInitialize:
    """``qclib``-style initialize gate for SPDMM state preparation.

    Mirrors the shape of ``qclib.state_preparation.LowRankInitialize`` so
    benchmark harnesses can swap SPDMM in seamlessly.

    Attributes
    ----------
    definition : QuantumCircuit
        The state-preparation circuit (output matches ``state_vector``
        up to a global phase).
    num_qubits : int
        Number of qubits, ``log2(len(state_vector))``.
    """

    def __init__(self, state_vector, label="spdmm"):
        v = np.asarray(state_vector, dtype=complex).reshape(-1)
        n = int(round(np.log2(len(v))))
        v = v / np.linalg.norm(v)
        self.definition = state_preparation(v)
        self.num_qubits = n
        self.label = label


__all__ = [
    "state_preparation",
    "siable",
    "siable_low_rank",
    "SiableBlockEncoding",
    "SpdmmInitialize",
    "SpdmmInitializeRust",
    "state_preparation_rust",
]

__version__ = "0.2.0"
