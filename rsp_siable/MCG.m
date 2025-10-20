function circuit = MCG( target, U, controls, controlStates, circuit )
% MCG() -- Simulate the qauntum circuit to prepare multiple controlled
% special unitary single-qubit gate in a non-strict way
% 
% You can refer --
% Efficient implementation of multicontrolled quantum gates, Phys. Rev. Appl., vol. 24, p. 044 030, 4 2025
%  or
%  Circuit decomposition of multicontrolled special unitary single-qubit gates. Trans. Comp.-Aided Des. Integ. Cir. Sys., 43(3):802–811, March 2024. 
% to get the optimized circuit with fewer CNOT gates. 
% ---------------------------------------- 
% input:
% target -- the target qubit 
% U -- the target single-qubit unitary 
% controls --  the control qubits 
% controlStates -- the control state in that control qubits 
% circuit -- the input circuit 
% ---------------------------------------- 
% ouput:
% circuit -- the output circuit 


    % ---------------------------------------------------------------------
    % Step 1. Compute parameter 
    % such that 
    % U = exp(1i*alpha) * [exp(-1i*beta/2), 0; 0, exp(1i*beta/2)] * [cos(gamma/2), -sin(gamma/2); sin(gamma/2), cos(gamma/2)] * [exp(-1i*delta/2), 0; 0, exp(1i*delta/2)];
    alpha = angle(det(U)) / 2 ; 
    if norm(U(1,1)) < 1e-14 
        beta = angle( -U(2,1) / U(1,2) ) ./ 2 ; 
        delta = -angle( -U(2,1) / U(1,2) ) ./ 2 ; 
    elseif norm(U(1,2)) < 1e-14 
        beta = angle( U(2,2) / U(1,1) ) ./ 2 ; 
        delta = angle( U(2,2) / U(1,1) ) ./ 2 ; 
    else 
        beta = ( angle( -U(2,1) / U(1,2) ) + angle( U(2,2) / U(1,1) )) ./ 2 ; 
        delta = ( -angle( -U(2,1) / U(1,2) ) + angle( U(2,2) / U(1,1) )) ./ 2 ; 
    end 
    gamma = angle(exp(1i*(-alpha+beta/2+delta/2)).*U(1,1) + 1j.*exp(1i*(-alpha-beta/2+delta/2)).*U(2,1)) .* 2 ; 
    
    % Step 2. U = exp(1j*alpha)*A*X*B*X*C, and eye(2) = A*B*C.
    % C = Rz((delta-beta)/2) ; 
    circuit.push_back( qclab.qgates.RotationZ(target, (delta-beta)/2) ) ; 
    % MCX 
    circuit.push_back( qclab.qgates.MCX(controls, target, controlStates) ) ; 
    % B = Ry(-gamma/2) * Rz(-(delta+beta)/2) ; 
    circuit.push_back( qclab.qgates.RotationZ(target, -(delta+beta)/2) ) ; 
    circuit.push_back( qclab.qgates.RotationY(target, -gamma/2) ) ; 
    % MCX 
    circuit.push_back( qclab.qgates.MCX(controls, target, controlStates) ) ; 
    % A = Rz(beta) * Ry(gamma/2) ; 
    circuit.push_back( qclab.qgates.RotationY(target, gamma/2) ) ; 
    circuit.push_back( qclab.qgates.RotationZ(target, beta) ) ;
    % exp(1j*alpha) 
    %% If 
    % circuit = MCPhase1( n, target, alpha, controls, controlStates, circuit ) ; 
end

