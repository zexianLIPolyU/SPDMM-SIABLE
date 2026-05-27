"""
Example 3 — Low-rank block encoding with automatic method dispatch.

Sweeps (n, K) pairs and reports which of the four constructions
{state_prep, isometry, spdmm_half, full_rank} was selected, together
with the CNOT count.  Reproduces the first three rows of Table V of
the paper.
"""

import numpy as np
from qiskit.quantum_info import Operator
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from spdmm import siable_low_rank


def _rank_K_matrix(N, K, rng):
    W = rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))
    V = rng.standard_normal((K, N)) + 1j * rng.standard_normal((K, N))
    A = W @ V
    A /= max(np.linalg.norm(A, 2), 1.0)         # ||A||_2 <= 1
    return A


def main():
    rng = np.random.default_rng(1)
    print(f"{'matrix N':>9} {'K':>3}   {'method':>12}   "
          f"{'CNOTs':>6}   {'encodes':>7}   {'block err':>10}")
    print("-" * 72)
    for n in [3, 4, 5, 6]:
        N = 2 ** n
        for K in [1, 2, 3]:
            if K >= N:
                continue
            A = _rank_K_matrix(N, K, rng)
            qc, alpha, info = siable_low_rank(A, K, return_info=True)
            U = Operator(qc).data
            block = alpha * U[0::2, 0::2]
            # Compare against A or its rank-K truncation, per info['encodes']
            target = A
            if info["encodes"] == "A_K":
                Uu, Sv, Vh = np.linalg.svd(A)
                target = Uu[:, :K] @ np.diag(Sv[:K]) @ Vh[:K, :]
            err = np.linalg.norm(target - block)
            print(f"{N:>9} {K:>3}   {info['method']:>12}   "
                  f"{info['nCNOT']:>6}   {info['encodes']:>7}   {err:>10.2e}")


if __name__ == "__main__":
    main()
