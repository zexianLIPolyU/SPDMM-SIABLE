function [alpha, gamma, beta, delta] = u2_decomposition(U) 
    % U2_DECOMPOSITION_ROBUST Decomposes a U(2) matrix into Rz(psi) * Ry(theta) * Rz(phi) with global phase
    % Input: U - 2x2 unitary matrix
    % Outputs: alpha - global phase
    %          beta, delta, gamma - rotation angles for Ry and Rz gates
    % such that 
    % ---------------------------------------------------------------------
    % circuit = qclab.QCircuit(1)
    % circuit.push_back( qclab.qgates.U3( 0, gamma, beta, delta ) )
    % norm( exp(1j*alpha) .* circuit.matrix - U ) 
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
    
    % Step 2. QCLAB U3 function global phase shift
    alpha = alpha - beta / 2 - delta / 2 ;

end