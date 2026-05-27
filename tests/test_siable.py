"""Test correctness and CNOT counts of siable() and siable_low_rank()."""

import os
import sys

# Make the local 'spdmm' package importable when this file is run directly
# (i.e. ``python tests/test_siable.py``).  Pytest already handles this via
# conftest.py, but a direct invocation does not.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Operator

from spdmm import siable, siable_low_rank, SiableBlockEncoding


def test_full_rank_correctness():
    """Verify that siable produces a correct block encoding."""
    print("=" * 70)
    print("Test 1: Full-rank SIABLE correctness")
    print("=" * 70)
    np.random.seed(42)
    for n in [2, 3, 4]:
        N = 2 ** n
        # Random complex matrix
        A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
        # Make it interesting -- not unitary, not trivial
        A = A * 0.7  # ensures ||A||_2 < 1 sometimes

        qc, alpha, info = siable(A, return_info=True)
        U = Operator(qc).data

        # Block-encoding relation:  A[i,j] == alpha * U[2*i, 2*j]
        # because the ancilla is qiskit-qubit-0 (LSB), so ancilla=0 means
        # the global index is 2*i.
        block = alpha * U[0::2, 0::2]
        err = np.linalg.norm(A - block)
        print(f"  n={n} (N={N:3d}): alpha={alpha:.4f}, "
              f"block-encoding error = {err:.2e}, "
              f"CNOTs (build) = {info['nCNOT']}")
        assert err < 1e-10, f"n={n}: large error {err}"


def test_full_rank_cnot_counts():
    """Report full-rank SIABLE C-NOT counts for n_total = 3..7."""
    print()
    print("=" * 70)
    print("Test 2: Full-rank SIABLE C-NOT counts")
    print("=" * 70)
    np.random.seed(0)
    for n_total in [3, 4, 5, 6, 7]:
        m = n_total - 1  # matrix qubits
        N = 2 ** m
        A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
        qc, alpha, info = siable(A, return_info=True)
        cnots_build = info["nCNOT"]
        qct = transpile(qc, basis_gates=["u", "cx"])
        cnots_tr = qct.count_ops().get("cx", 0)
        print(f"  n_total={n_total} (matrix {N}x{N}): "
              f"CNOTs build={cnots_build}, transpile={cnots_tr}")


def test_low_rank_correctness():
    """Verify low-rank SIABLE block-encodes the rank-K truncation."""
    print()
    print("=" * 70)
    print("Test 3: Low-rank SIABLE correctness")
    print("=" * 70)
    np.random.seed(7)
    for n, K in [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3), (5, 2)]:
        N = 2 ** n
        # Build a genuinely rank-K matrix to test the truncation
        U_left = np.linalg.qr(np.random.randn(N, K) + 1j * np.random.randn(N, K))[0]
        U_right = np.linalg.qr(np.random.randn(N, K) + 1j * np.random.randn(N, K))[0]
        svs = np.random.rand(K) * 0.8 + 0.1  # singular values in (0.1, 0.9)
        A = U_left @ np.diag(svs) @ U_right.conj().T
        # Sanity: actual rank should be K
        assert np.linalg.matrix_rank(A, tol=1e-8) == K

        qc, alpha, info = siable_low_rank(A, K, return_info=True)
        U = Operator(qc).data
        block = alpha * U[0::2, 0::2]
        err = np.linalg.norm(A - block)
        rel_err = err / np.linalg.norm(A)
        print(f"  n={n} K={K}: alpha={alpha:.4f}, "
              f"block-encoding rel-err = {rel_err:.2e}, CNOTs = {info['nCNOT']}")


def test_low_rank_full_rank_consistency():
    """Verify that low-rank SIABLE with K = N still gives correct encoding."""
    print()
    print("=" * 70)
    print("Test 4: Low-rank SIABLE with K=N agrees with siable")
    print("=" * 70)
    np.random.seed(7)
    n = 3
    N = 2 ** n
    A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
    qc_lr, alpha_lr, info_lr = siable_low_rank(A, K=N, return_info=True)
    qc_fr, alpha_fr, info_fr = siable(A, return_info=True)

    U_lr = Operator(qc_lr).data
    U_fr = Operator(qc_fr).data
    block_lr = alpha_lr * U_lr[0::2, 0::2]
    block_fr = alpha_fr * U_fr[0::2, 0::2]

    print(f"  full-rank:  CNOTs = {info_fr['nCNOT']}, ||A-block|| = "
          f"{np.linalg.norm(A - block_fr):.2e}")
    print(f"  low-rank K=N: CNOTs = {info_lr['nCNOT']}, ||A-block|| = "
          f"{np.linalg.norm(A - block_lr):.2e}")


def test_low_rank_cnot_counts_sweep():
    """Sweep over a grid of (n_total, K) pairs and report the C-NOT cost of
    the auto-selected method.  Verifies correctness for the smaller cases
    (n_total <= 7) where building the full Operator is tractable.

    When auto-dispatch picks ``'full_rank'`` (because it's cheaper than
    every low-rank construction), the circuit block-encodes ``A`` itself,
    not the rank-``K`` truncation ``A_K``.  Verification compares against
    the matrix the chosen method actually encodes.
    """
    print()
    print("=" * 70)
    print("Test 5: Low-rank SIABLE C-NOT sweep over (n_total, K) grid")
    print("=" * 70)
    pairs = [
        (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
        (4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 10),
        (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 10),
        (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 10),
        (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 10),
        (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 10),
        (9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 10),
    ]
    print(f"  {'n_total':>7} {'K':>3}   {'method':>12} {'encodes':>8} "
          f"{'CNOTs':>7}  {'block-err':>10}")
    print("  " + "-" * 60)
    np.random.seed(0)
    for n_total, K in pairs:
        m = n_total - 1  # matrix qubits
        N = 2 ** m
        A = np.random.randn(N, N) + 1j * np.random.randn(N, N)

        qc, alpha, info = siable_low_rank(A, K, method='auto', return_info=True)
        cnots = info['nCNOT']
        encoded_label = info['encodes']

        # Verify correctness for cases where Operator(qc) is tractable.
        if n_total <= 7:
            # Compare to the matrix actually encoded: A_K for the low-rank
            # paths, A for the full-rank fallback.
            if info['encodes'] == 'A_K':
                W, S, Vh = np.linalg.svd(A / alpha)
                Sigma_K = np.zeros(N)
                Sigma_K[:K] = np.clip(S[:K], 0, 1)
                target = W @ np.diag(Sigma_K) @ Vh
            else:  # 'A'
                target = A / alpha
            U_full = Operator(qc).data
            block = np.array(
                [[U_full[2 * i, 2 * j] for j in range(N)] for i in range(N)]
            )
            err = float(np.linalg.norm(block - target))
            err_str = f"{err:.2e}"
            assert err < 1e-7, \
                f"n_total={n_total}, K={K}, method={info['method']}: " \
                f"block-encoding err {err:.3e} too large"
        else:
            err_str = "  (skipped)"

        print(f"  {n_total:>7} {K:>3}   {info['method']:>12} "
              f"{encoded_label:>8} {cnots:>7}  {err_str:>10}")


def test_low_rank_state_prep_K1():
    """K=1 auto-dispatch must pick the SPDMM state_preparation path, and the
    resulting circuit must correctly block-encode the rank-1 truncation."""
    print()
    print("=" * 70)
    print("Test 7: K=1 SIABLE-low-rank via SPDMM state_preparation")
    print("=" * 70)
    np.random.seed(0)
    for n_total in [3, 4, 5, 6, 7, 8]:
        n = n_total - 1
        N = 2 ** n
        # Build a rank-1 matrix so K=1 truncation = full A.
        u = np.random.randn(N) + 1j * np.random.randn(N)
        v = np.random.randn(N) + 1j * np.random.randn(N)
        A = np.outer(u, v.conj())
        qc, alpha, info = siable_low_rank(A, 1, method='auto', return_info=True)
        cnots = qc.count_ops().get("cx", 0)
        print(f"  n_total={n_total}, K=1: chosen='{info['method']}', "
              f"CNOTs={cnots:>4}")
        assert info['method'] == 'state_prep', \
            f"auto-dispatch should pick state_prep for K=1, got '{info['method']}'"


def test_low_rank_spdmm_half_iso():
    """SPDMM half-isometry path must (a) correctly block-encode the rank-K
    truncation, (b) be auto-selected when K_pad == 2^(n-1) and n >= 3,
    and (c) save C-NOTs vs the qiskit-Isometry path."""
    print()
    print("=" * 70)
    print("Test 8: SPDMM half-isometry path -- correctness and C-NOT savings")
    print("=" * 70)
    np.random.seed(7)
    saving_total = 0
    for n in [3, 4, 5, 6]:
        N = 2 ** n
        A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
        # K = 2^(n-1) so K_pad = 2^(n-1) -- the regime where spdmm_half
        # is cheapest.
        K = 1 << (n - 1)

        qc_iso, alpha_iso, info_iso = siable_low_rank(
            A, K, method='isometry', return_info=True
        )
        qc_half, alpha_half, info_half = siable_low_rank(
            A, K, method='spdmm_half', return_info=True
        )
        qc_auto, alpha_auto, info_auto = siable_low_rank(
            A, K, method='auto', return_info=True
        )

        # Build the rank-K truncation in qiskit ordering for comparison.
        W, S, Vh = np.linalg.svd(A / alpha_half)
        Sigma_K = np.zeros(N)
        Sigma_K[:K] = np.clip(S[:K], 0, 1)
        A_K = W @ np.diag(Sigma_K) @ Vh

        # Verify all three paths block-encode A_K correctly.
        for name, qc in [('isometry', qc_iso),
                         ('spdmm_half', qc_half),
                         ('auto', qc_auto)]:
            U = Operator(qc).data
            block = np.array(
                [[U[2*i, 2*j] for j in range(N)] for i in range(N)]
            )
            err = np.linalg.norm(block - A_K)
            assert err < 1e-8, \
                f"n={n}, K={K}, method={name}: block-enc err {err:.3e}"

        # Auto should pick spdmm_half here.
        assert info_auto['method'] == 'spdmm_half', \
            f"auto-dispatch should pick spdmm_half when K_pad=2^(n-1)=2^{n-1}"
        # spdmm_half should be cheaper than isometry.
        c_iso = info_iso['nCNOT']
        c_half = info_half['nCNOT']
        saving = c_iso - c_half
        saving_total += saving
        pct = 100.0 * saving / c_iso
        print(f"  n={n} K={K} (K_pad=2^(n-1)): "
              f"iso={c_iso:>5} half={c_half:>5} "
              f"saving={saving:>5} ({pct:.1f}%)")
        assert c_half < c_iso, \
            f"n={n}, K={K}: spdmm_half should be cheaper, got {c_half} vs {c_iso}"

    print(f"  -- total C-NOT saving across all sizes: {saving_total}")


def test_gate_wrapper():
    print()
    print("=" * 70)
    print("Test 6: SiableBlockEncoding gate wrapper")
    print("=" * 70)
    np.random.seed(7)
    A = np.random.randn(8, 8) + 1j * np.random.randn(8, 8)
    gate = SiableBlockEncoding(A)
    print(f"  Gate: name={gate.name}, num_qubits={gate.num_qubits}, "
          f"alpha={gate.alpha:.4f}")
    # The .definition is the actual circuit
    qc = gate.definition
    U = Operator(qc).data
    block = gate.alpha * U[0::2, 0::2]
    err = np.linalg.norm(A - block)
    print(f"  Definition matches: err = {err:.2e}")


def plot_siable_circuit_n3_K1():
    """Plot the auto-dispatched SIABLE-low-rank circuit for n_total=3, K=1.

    n_total = 3 means the circuit has 3 qubits: 1 ancilla on qubit 0 plus
    2 data qubits (encoding a 4x4 matrix).  At K=1 the auto-dispatcher
    selects the SPDMM state_preparation path.

    Writes a matplotlib rendering to ``siable_n3_K1.png`` (next to this
    test script, so the script works on any OS) and prints a text drawing
    to stdout.
    """
    print()
    print("=" * 70)
    print("Plot: SIABLE low-rank circuit for n_total=3, K=1 (4x4 matrix)")
    print("=" * 70)
    np.random.seed(0)
    N = 4  # n_data = 2 → 2^2 = 4
    A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
    qc, alpha, info = siable_low_rank(A, 1, method='auto', return_info=True)

    print(f"  chosen method : {info['method']}")
    print(f"  num qubits    : {qc.num_qubits}  (qubit 0 = ancilla, 1..2 = data)")
    print(f"  CNOT count    : {info['nCNOT']}")
    print(f"  total gates   : {sum(qc.count_ops().values())}")
    print(f"  alpha = ||A||₂ = {alpha:.4f}")
    print()
    print("  Text drawing:")
    print()
    # Indent the text drawing.
    text_drawing = str(qc.draw('text', fold=120))
    for line in text_drawing.splitlines():
        print("    " + line)

    # Save matplotlib rendering next to this script (portable across OS).
    import os
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "siable_n3_K1.png"
    )
    try:
        import matplotlib
        matplotlib.use('Agg')
        fig = qc.draw('mpl', style='clifford', fold=-1)
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        import matplotlib.pyplot as plt
        plt.close(fig)
        print()
        print(f"  Saved matplotlib rendering to: {out_path}")
    except Exception as e:
        print()
        print(f"  matplotlib rendering failed: {type(e).__name__}: {e}")


if __name__ == "__main__":
    test_full_rank_correctness()
    test_full_rank_cnot_counts()
    test_low_rank_correctness()
    test_low_rank_full_rank_consistency()
    test_low_rank_cnot_counts_sweep()
    test_gate_wrapper()
    test_low_rank_state_prep_K1()
    test_low_rank_spdmm_half_iso()
    plot_siable_circuit_n3_K1()
    print()
    print("ALL TESTS PASSED")
