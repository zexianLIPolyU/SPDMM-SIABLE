"""
Example 1 — n-qubit state preparation via SPDMM.

Builds a state-preparation circuit that maps |0...0> to a target state psi,
verifies the output state matches psi, and reports the C-NOT count.
"""

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Statevector

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from spdmm import state_preparation


def main():
    rng = np.random.default_rng(42)
    for n in [3, 4, 5, 6]:
        N = 2 ** n
        psi = rng.standard_normal(N) + 1j * rng.standard_normal(N)
        psi /= np.linalg.norm(psi)

        qc = state_preparation(psi)
        prepared = Statevector.from_instruction(qc).data
        # SPDMM matches up to a global phase
        overlap = np.vdot(prepared, psi)
        phase = np.conj(overlap) / abs(overlap)
        err = np.linalg.norm(prepared - phase * psi)

        qc_compiled = transpile(qc, basis_gates=["u", "cx"], optimization_level=0)
        n_cx = qc_compiled.count_ops().get("cx", 0)

        print(f"n={n}: ||psi_out - psi|| = {err:.2e}, CNOTs = {n_cx}")


if __name__ == "__main__":
    main()
