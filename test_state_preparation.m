% test state_preparation
clc ; clear ; close all ; 
addpath('QCLAB');
addpath("rsp_siable");


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
