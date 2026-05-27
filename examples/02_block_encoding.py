"""
Example 2 — Full-rank single-ancilla block encoding via SIABLE.

Block-encodes a random 2^n x 2^n matrix with spectral norm <= 1 and
verifies (i) the block-encoding relation A == alpha * <0_anc|U|0_anc>
and (ii) the CNOT count matches the paper Table III.
"""

import numpy as np
from qiskit.quantum_info import Operator
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from spdmm import siable


def main():
    # Paper Table III, n = total qubits (matrix is 2^{n-1} x 2^{n-1}).
    paper_cnots = {3: 9, 4: 45, 5: 205, 6: 877, 7: 3629}

    rng = np.random.default_rng(0)
    for n_paper in [3, 4, 5]:
        N = 2 ** (n_paper - 1)                  # matrix is N x N = 2^{n-1} x 2^{n-1}
        A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        A *= 0.3                                # ||A||_2 not too close to 1

        qc, alpha, info = siable(A, return_info=True)
        assert qc.num_qubits == n_paper         # 1 ancilla + (n_paper - 1) data

        U = Operator(qc).data
        block = alpha * U[0::2, 0::2]
        err = np.linalg.norm(A - block)

        match = "OK" if info["nCNOT"] == paper_cnots[n_paper] else "(differs)"
        print(f"n={n_paper} (matrix {N}x{N}): ||A - alpha*block|| = {err:.2e}, "
              f"CNOTs = {info['nCNOT']}  (paper: {paper_cnots[n_paper]}) {match}")


if __name__ == "__main__":
    main()
