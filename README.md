# Recursive State Preparation (RSP) and Single Ancilla Block Encoding (SIABLE)
We present a MATLAB implementation of the Recursive State Preparation (RSP) and Single Ancilla Block Encoding (SIABLE) protocols, designed to use as few CNOT gates as possible. This implementation is built upon the [QCLAB](https://github.com/QuantumComputingLab/qclabs) framework.

## Problem Formulation
1.  **State Preparation:** Given $n \in Z_+$ and a vector $$\[\psi_i\]_{i=0}^{2^n-1}$$ whose $$\ell_2$$-norm is $$1$$, generate the quantum state:

$$
    \ket{\psi} = \sum_{i=0}^{2^n-1} \psi_i |i\rangle
$$
    
from the initial state $|0\rangle^{\otimes n}$.
    
2. **Block Encoding:** Given a complex matrix $A$, generate a quantum circuit whose left-upper block of the matrix form $U$ is $A/\alpha$ as
   
$$
    U = \begin{bmatrix} A/\alpha & * \\
    * & *  \end{bmatrix}.
$$

Formally, a $(\alpha,a,\epsilon)$-block-encoding statisfy that 

$$
   \Vert (\bra{0}^{\otimes a}\otimes I)  U  (\ket{0}^{\otimes a}\otimes I) - A/\alpha \Vert \leq \varepsilon,
$$

where $\alpha$ is the subnotmalization, $a$ is the number of ancillas, $\varepsilon$ is the precision. The optimal subnotmalization of a block encoding is the spectral norm $\Vert A\Vert_2$ and the fewest number of ancillas is $1$. 


**Question:**

**How to use as few CNOT gates as possible to generate quantum circuit of state preparation and ($$\Vert A\Vert_2$$,1,0)-block-encoding ?**


# Implementation 

In order to run the MATLAB implementation of Recursive State Preparation (RSP) and Single Ancilla Block Encoding Protocol (SIABLE):
1. Down [rsp-siable](https://github.com/zexianLIPolyU/RSP-SIABLE/tree/main/rsp_siable) and [QCLAB](https://github.com/QuantumComputingLab/qclab) repositories.
2. Unzip it and add `rsp_siable` and `QCLAB` files into your MATLAB path.
    ```
    cd("rsp_siable")
    cd("QCLAB")
    ```
3. State preparation can be run by
     ```
      n = 5;
      N = pow2(n) ; 
      state_complex = randn(N,1) + randn(N,1) .* 1j ; 
      state_complex = state_complex ./ norm(state_complex,2) ; 
      logging = true; % no record 
      [circuit, global_phase, CNOT_count] = state_preparation( state_complex, 1, logging ) ; 
      fprintf("N_state(%d) = %d\n\n ", n, CNOT_count) ; 
      circuit.draw()
      M = circuit.matrix;
      norm(M(:,1) .* global_phase - state_complex )
     ```

4. Test single ancilla block encoding protocol for full-rank and low-rank matrix 
   
      See demo in [test_siable_image.mlx](https://github.com/zexianLIPolyU/RSP-SIABLE/blob/main/test_siable_image.mlx). 






# Result 
## 1. Recursive State Preparation (RSP)

**Table: Comparison of the number of C-NOT gates between proposed recursive state preparation method (RSP) and other state preparation algorithms.**

| Methods | Script | 2 | 3 | 4 | 5 | 10 | 15 | Leading constant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PB | [qclib](https://github.com/qclib/qclib) | 1 | 4 | 9 | 26 | 919 | 38813 | $23/24$ |
| Isometry | [qiskit](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.StatePreparation) | 1 | 4 | 11 | 26 | 1013 | 32752 | $23/24$ |
| LRSP | [qclib](https://github.com/qclib/qclib) | 1 | 4 | 9 | 21 | 913 | 30999 | $23/24$ |
| **RSP** **(Proposed method)** | [rsp](https://github.com/zexianLIPolyU/siable/blob/main/test_state_preparation.m) | 1 | **3** | **7** | **18** | **867** | **29627** | $11/12$ |
| **Theoretical lower bound for state preparation**| - | 1 | 2 | 5 | 12 | 505 | 16373 | $1/2$

## 2. Single Ancilla Block Encoding Protocol (SIABLE)

### 2.1. Full-rank matrix encoding 

Script: [siable](https://github.com/zexianLIPolyU/RSP-SIABLE/blob/main/rsp_siable/siable.m)

**Table: Comparison of the number of C-NOT gates between the single ancilla block encoding protocol (SIABLE) for general $2^{n-1}\times 2^{n-1}$ full-rank matrix and other unitary synthesis protocol and bounds in an $n$-qubit system.**

| Number of qubits | Script | 3 | 4 | 5 | 6 | 7 | n |
| --- | --- | --- | --- | --- | --- | --- | --- |
| QSD | [qiskit](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.transpiler.passes.UnitarySynthesis) | 20 | 100 | 444 | 1868 | 7660 | $(23/48)\times4^n - (3/2)\times 2^n + (4/3)$ |
| Block-ZXZ | [siable](https://github.com/zexianLIPolyU/siable/blob/main/test_siable_CNOT.m) | 19 | 95 | 423 | 1783 | 7319 | $(22/48)\times4^n - (3/2)\times 2^n + (5/3)$ |
| Shende's lower bound | - | 14 | 61 | 252 | 1020 | 4091 | $\lceil (1/4)\times(4^n - 3n - 1) \rceil$ |
| **SIABLE for full-rank matrix** <br> **(Proposed method)** | [siable](https://github.com/zexianLIPolyU/siable/blob/main/test_siable_CNOT.m) | **9** | **45** | **205** | **877** | **3629** | $(11/48)\times 4^n - 2^n + (7/3)$ |
| **Theoretical lower bound for block encoding**| - | 6 | 29 | 125 | 508 | 2043 | $\lceil (1/8)\times4^n - (3/4)\times n \rceil$ 

### 2.2 Low-rank matrix encoding

Script: [siable_low_rank](https://github.com/zexianLIPolyU/RSP-SIABLE/blob/main/rsp_siable/siable_low_rank.m)

**Table: Number of C-NOT gates in the single ancilla block encoding protocol (SIABLE) for general low-rank and full-rank $2^{n-1}\times 2^{n-1}$ matrix with a single ancilla.**

| n\rank | 1   | 2    | 3    | 4    | 5    | 10    | full-rank |
|--------------------|-----|------|------|------|------|-------|-----------|
| 3                  | 6   |      |      |      |      |       | 9         |
| 4                  | 14  |      |      |      |      |       | 45        |
| 5                  | 30  |      |      |      |      |       | 205       |
| 6                  | 68  | 616  |      |      |      |       | 877       |
| 7                  | 148 | 662  | 1098 | 1532 | 1970 |       | 3629      |
| 8                  | 314 | 1178 | 1940 | 2700 | 3464 | 7064  | 14765     |
| 9                  | 654 | 2044 | 3308 | 4570 | 5836 | 11898 | 59565


