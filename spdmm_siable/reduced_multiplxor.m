function [ circuit, global_phase, R, nCNOT ] = reduced_multiplxor( n, U2_sequence, controls, target, circuit ) 
% REDUCED_MULTIPLXOR --
% Decouple k-fold uniformly controlled U(2) gate with 2^k-1 CNOT gates.
% ------------------------------------------------------------------------
% Generate circuit to implement recursive decomposition of the k-fold 
% uniformly controlled U(2) gate, eliminating the uniform control node on the
% qubit m as 
% 
%       F_t^k[U(2)] = F_m^t[R_z] F_t^{k-1}[U(2)] D_t^m F_t^{k-1}[U(2)] 
% 
% reference: Quantum circuits with uniformly controlled one-qubit gates
% ------------------------------------------------------------------------
% test program:
% n = 5 ;
% U2_sequence = [] ; 
% for i = 1:pow2(n-1) 
%     U2_sequence = cat(3, U2_sequence, generate_U2n(1)) ; 
% end 
% [ circuit, global_phase, R ] = multiplxor( n, U2_sequence, 0:n-2, n-1) ;
% M = global_phase .* diag(R) * circuit.matrix ; 
% err = 0 ; 
% for i = 1:pow2(n-1) 
%     err = err + norm(M(2*i-1:2*i,2*i-1:2*i) - U2_sequence(:,:,i)) ; 
% end
% err 
% ------------------------------------------------------------------------

    if nargin <= 4 
        circuit = qclab.QCircuit( target + 1 ) ; 
    end
    assert( n >= 2 ) ; 
    % pow2(control nodes) is equal to the size of U2_sequence
    assert( size(U2_sequence,3) == pow2(size(controls,2)) ) ; 
    if size(controls,2) == 0 
    % Independent case, implement a single-qubit gate
        [alpha, gamma, beta, delta] = u2_decomposition(U2_sequence) ; 
        circuit.push_back( qclab.qgates.U3(target, gamma, beta, delta ) ) ; 
        global_phase = exp(1j.*alpha) ; 
        R = false ; 
        return ; 
    end 

    A_U2_sequence = [] ; 
    B_U2_sequence = [] ; 
    N = pow2(size(controls,2)) ; 
    R = zeros( 1, 2*N ) ; 
    nCNOT = 0 ; 
    
    for i = 1 : N/2 
        a = U2_sequence(:,:,i) ; 
        b = U2_sequence(:,:,i + N/2) ; 
        X_ab = a * b' ; 
        delta = pi / 2 ; 
        phi = angle( det(X_ab) ) ; 
        x1 = X_ab(1,1) .* exp(-1j.*phi./2) ; 
        r1 = exp( 1j./2 .* (delta - phi./2 - angle(x1)) ) ; 
        r2 = exp( 1j./2 .* (delta - phi./2 + angle(x1) + pi) ) ; 
        r = diag([r1, r2]) ; 
        X_ab = r * X_ab * r ; 
        [u, d2] = eig(X_ab) ; 
        d = sqrt(d2) ; 
        v = d * u' * r' * b ; 
        A_U2_sequence = cat(3, A_U2_sequence, v) ; 
        B_U2_sequence = cat(3, B_U2_sequence, u) ; 
        R(2*i-1:2*i) = [r1',r2'] ; 
        R(2*i-1+N:2*i+N) = [r1,r2] ; 
    end 
    if size(controls,2) == 1 
        [alpha0, gamma, beta, delta] = u2_decomposition( A_U2_sequence ) ; 
        circuit.push_back( qclab.qgates.U3( target, gamma, beta, delta ) ) ; 
        [circuit_D, global_phase_D] = circuit_diag( controls(1), target ) ; 
        nCNOT = nCNOT + 1 ; 
        circuit.push_back( circuit_D ) ; 
        [alpha1, gamma, beta, delta] = u2_decomposition( B_U2_sequence ) ; 
        circuit.push_back( qclab.qgates.U3( target, gamma, beta, delta ) ) ; 
        global_phase = exp(1i.*(alpha0+alpha1)) * global_phase_D ; 
    else
        [ circuitA, global_phase0, RA, nCNOT0 ] = reduced_multiplxor( n, A_U2_sequence, controls(2:end), target ) ; 
        nCNOT = nCNOT + nCNOT0 ; 
        circuit.push_back( circuitA ) ; 
        [ circuit_D, global_phase_D ] =  circuit_diag( controls(1), target ) ; 
        nCNOT = nCNOT + 1 ; 
        circuit.push_back( circuit_D ) ; 
        for i = 1:N/2 
            B_U2_sequence(:, :, i) = B_U2_sequence(:, :, i) .* RA(2*i-1:2*i) ; 
        end 
        [ circuitB, global_phase1, RB, nCNOT0 ] = reduced_multiplxor( n, B_U2_sequence, controls(2:end), target ) ; 
        nCNOT = nCNOT + nCNOT0 ; 
        circuit.push_back( circuitB ) ;  
        R = R .* kron([1,1], RB) ; 
        global_phase = global_phase0 * global_phase_D * global_phase1 ; 
    end
end

function [circuit, global_phase] = circuit_diag(control, target)
% generate a quantum circuit with diagonal matrix
% expm(1j.*pi./4.*kron(diag([1,-1]),diag([1,-1]))) 
% ------------------------------------------------
% test program: 
    % [circuit, global_phase] = circuit_diag(0, 1) ;
    % norm( circuit.matrix.* global_phase - expm(1j.*pi./4.*kron(diag([1,-1]),diag([1,-1])))  ) 
% ------------------------------------------------
    circuit = qclab.QCircuit( target + 1 ) ;
    circuit.push_back( qclab.qgates.Hadamard(target) ) ; 
    circuit.push_back( qclab.qgates.CNOT(control, target) ) ; 
    circuit.push_back( qclab.qgates.Hadamard(target) ) ; 
    circuit.push_back( qclab.qgates.RotationZ(control, -pi/2) ) ; 
    circuit.push_back( qclab.qgates.RotationZ(target, -pi/2) ) ; 
    global_phase = exp(-1i.*pi/4) ; 
end