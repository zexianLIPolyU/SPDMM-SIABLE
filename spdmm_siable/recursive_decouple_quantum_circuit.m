function [circuit, global_phase, Delta, CNOT_count] = recursive_decouple_quantum_circuit(circuit, U, n, global_phase, Delta, is_last_U, logging ) 
% Recursive decouple a pow2(n)*pow2(n) unitary U 
% Input: 
%   circuit      -- qclab circuit 
%   U            -- unitary to be decoupled 
%   n            -- qubit number to represent U 
%   global_phase -- the global phase of the circuit 
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict 
%   is_last_U    -- boolean logical that whether U is the last unitary in the decomposition 
% --------------------------------------------------- 
% Output: 
%   circuit      -- qclab circuit 
%   global_phase -- the global phase of the circuit 
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict 
% ---------------------------------------------------
% test program: 
    % circuit = qclab.QCircuit( n );
    % [circuit, global_phase, Delta, CNOT_count] = recursive_decouple_quantum_circuit( circuit, U, n, 1, eye(4), true ) ; 
    % norm( circuit.matrix .* global_phase - U ) 
% ---------------------------------------------------
    
    assert( (pow2(n) == size(U, 1)) & (pow2(n) == size(U, 2)) )
    CNOT_count = 0 ; 
    if n > 2 
        % Step 1. Compute the rotation parameters 
        [A, B1, B2, C] = transfomed_ZXZ( U ) ; 
        % norm( U - [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), A] * kron(H, eye(pow2(n-1))) * [B1, zeros(pow2(n-1)); zeros(pow2(n-1)), B2] * kron(H, eye(pow2(n-1))) * [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), C] )
        [UA, DA, VA] = decouple_unitary( eye(pow2(n-1)), A ) ; 
        % norm([UA, zeros(pow2(n-1)); zeros(pow2(n-1)), UA] * [DA, zeros(pow2(n-1)); zeros(pow2(n-1)), DA'] * [VA, zeros(pow2(n-1)); zeros(pow2(n-1)), VA] - [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), A] ) 
        [UC, DC, VC] = decouple_unitary( eye(pow2(n-1)), C ) ; 
        % norm([UC, zeros(pow2(n-1)); zeros(pow2(n-1)), UC] * [DC, zeros(pow2(n-1)); zeros(pow2(n-1)), DC'] * [VC, zeros(pow2(n-1)); zeros(pow2(n-1)), VC] - [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), C] ) 
        Z = [1 0; 0 -1] ; 
        [UB, DB, VB] = decouple_unitary( VA*B1*UC, kron(Z,eye(pow2(n-2)))*VA*B2*UC*kron(Z,eye(pow2(n-2))) ) ; 
        % norm([UB, zeros(pow2(n-1)); zeros(pow2(n-1)), UB] * [DB, zeros(pow2(n-1)); zeros(pow2(n-1)), DB'] * [VB, zeros(pow2(n-1)); zeros(pow2(n-1)), VB] - [VA*B1*UC, zeros(pow2(n-1)); zeros(pow2(n-1)), kron(Z,eye(pow2(n-2)))*VA*B2*UC*kron(Z,eye(pow2(n-2)))] )
        
        % Step 2. Construct quantum circuit 
        [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, VC, n-1, global_phase, Delta, false, logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        % M = kron(eye(2), VC) ;
        % norm( kron( eye(pow2(n-2)), Delta') * circuit.matrix .* global_phase -  M )  
        [circuit, CNOT_count0] = uniformly_controlled_rotation_z(circuit, DC, n, 'L', logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        % M = [DC, zeros(pow2(n-1)); zeros(pow2(n-1)), DC'] * M ; 
        % norm( kron( eye(pow2(n-2)), Delta') * kron(CNOT21, eye(pow2(n-2)))  * circuit.matrix .* global_phase -  M )  
        circuit.push_back( qclab.qgates.Hadamard(circuit.nbQubits - n) ) ; 
        [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, VB, n-1, global_phase, Delta, false, logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        [circuit, CNOT_count0] = uniformly_controlled_rotation_z(circuit, DB, n, 'M', logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, UB, n-1, global_phase, Delta, false, logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        circuit.push_back( qclab.qgates.Hadamard(circuit.nbQubits - n) ) ; 
        % M = kron(H, UB) * [DB, zeros(pow2(n-1)); zeros(pow2(n-1)), DB'] * kron(H, VB) * kron(CNOT21, eye(pow2(n-2))) * M ; 
        % norm( kron( eye(pow2(n-2)), Delta') * circuit.matrix .* global_phase -  M )  
        [circuit, CNOT_count0] = uniformly_controlled_rotation_z(circuit, DA, n, 'R', logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        % M = [DA, zeros(pow2(n-1)); zeros(pow2(n-1)), DA'] * kron(CNOT21, eye(pow2(n-2))) * M ; 
        % norm( kron( eye(pow2(n-2)), Delta') * circuit.matrix .* global_phase - M ) 
        [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, UA, n-1, global_phase, Delta, is_last_U, logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
        % M = kron(eye(2), UA) * M ; 
        % norm( kron( eye(pow2(n-2)), Delta') * circuit.matrix .* global_phase -  M )  
    else 
        [circuit, global_phase, Delta, CNOT_count0] = decoupleU4_quantum_circuit(circuit, U, global_phase, Delta, is_last_U, logging) ; 
        if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
    end 

end 



function [circuit, global_phase, Delta, CNOT_count] = decoupleU4_quantum_circuit(circuit, U, global_phase0, Delta0, is_last_U4, logging) 
% Decouple a 4*4 unitary U into a quantum circuit
% Input:
%   circuit      -- qclab circuit
%   U            -- unitary to be decoupled
%   global_phase -- the global phase of the circuit
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict
%   is_last_U4   -- boolean logical that whether U is the last U4 unitary in the decomposition
% ---------------------------------------------------
% Output:
%   circuit      -- qclab circuit 
%   global_phase -- the global phase of the circuit
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict
% ---------------------------------------------------
% test program: 
    % A = randn(4,4) + randn(4,4) .* 1j ; 
    % [A,~] = qr(A) ; 
    % circuit = qclab.QCircuit(2) ; 
    % global_phase0 = exp( 1j.*randn(1) ) ; 
    % Delta0 = expm( 1j.* diag(randn(4,1)) ) ; 
    % [circuit, global_phase, Delta] = decoupleU4_quantum_circuit( circuit, A, global_phase0, Delta0, false, false ) ; 
    % M = circuit.matrix ; 
    % norm( A * Delta0' .* global_phase0 - Delta' *  M(1:4,1:4) .* global_phase ) 
    % 
    % circuit = qclab.QCircuit(2) ; 
    % [circuit, global_phase] = decoupleU4_quantum_circuit( circuit, A, global_phase0, Delta0, true, false ) ; 
    % M = circuit.matrix ; 
    % norm( A * Delta0' .* global_phase0  - M(1:4,1:4) .* global_phase ) 
% ---------------------------------------------------
    
    if ~is_last_U4 

        [ LU1, LU2, RU1, RU2, thetaU, phiU, psiU, global_phaseU ] = KAK_decomposition2( U * Delta0' .* global_phase0 ) ; 
        [ alphaLU1, gammaLU1, betaLU1, deltaLU1 ] = u2_decomposition( LU1 ) ; 
        [ alphaLU2, gammaLU2, betaLU2, deltaLU2 ] = u2_decomposition( LU2 ) ; 
        [ alphaRU1, gammaRU1, betaRU1, deltaRU1 ] = u2_decomposition( RU1 ) ; 
        [ alphaRU2, gammaRU2, betaRU2, deltaRU2 ] = u2_decomposition( RU2 ) ; 
        
        circuit0 = qclab.QCircuit( 2, circuit.nbQubits - 2 ) ; 
        circuit0.push_back( qclab.qgates.U3(0, gammaLU1, betaLU1, deltaLU1 ) ) ; 
        circuit0.push_back( qclab.qgates.U3(1, gammaLU2, betaLU2, deltaLU2 ) ) ; 
        circuit0.push_back( qclab.qgates.CNOT(0,1) ) ; 
        circuit0.push_back( qclab.qgates.RotationX(0, thetaU) ) ; 
        circuit0.push_back( qclab.qgates.RotationZ(1, phiU) ) ; 
        circuit0.push_back( qclab.qgates.CNOT(0,1) ) ; 
        circuit0.push_back( qclab.qgates.U3(0, gammaRU1, betaRU1, deltaRU1 ) ) ; 
        circuit0.push_back( qclab.qgates.U3(1, gammaRU2, betaRU2, deltaRU2 ) ) ; 
        circuit.push_back(circuit0) ; 
        
        global_phase = 1 / global_phaseU * exp( 1j.*( alphaLU1 + alphaLU2 + alphaRU1 + alphaRU2 ) ) ; 
        Delta = diag( [exp(-1i.*psiU./2),exp(1i.*psiU./2),exp(1i.*psiU./2),exp(-1i.*psiU./2)] ) ; 

        if logging 
            CNOT_count = 2 ; 
        end 
        
        % M = circuit.matrix ;
        % norm( U * Delta0' .* global_phase0 - Delta' *  M(1:4,1:4) .* global_phase ) 

    else

        [ LU1, LU2, RU1, RU2, phiU, psiU, thetaU, global_phaseU ] = KAK_decomposition3( U * Delta0' .* global_phase0 ) ;
        [ alphaL1, gammaL1, betaL1, deltaL1 ] = u2_decomposition( LU1 ) ; 
        [ alphaL2, gammaL2, betaL2, deltaL2 ] = u2_decomposition( LU2 ) ; 
        [ alphaR1, gammaR1, betaR1, deltaR1 ] = u2_decomposition( RU1 ) ; 
        [ alphaR2, gammaR2, betaR2, deltaR2 ] = u2_decomposition( RU2 ) ; 
        global_phase = 1 / global_phaseU * exp( 1j * ( alphaL1 + alphaL2 + alphaR1 + alphaR2 )) ; 
        % global_phase_angle = -2 * angle( 1/ global_phaseU ) -2 .* ( alphaL1 + alphaL2 + alphaR1 + alphaR2 ) ;
        % global_phase = - 2 .* ( alphaLV1 + alphaLV2 + alphaRV1 + alphaRV2 ) ;

        circuit0 = qclab.QCircuit( 2, circuit.nbQubits - 2  ) ; 
        circuit0.push_back( qclab.qgates.U3(0, gammaL1, betaL1, deltaL1 ) ) ; 
        circuit0.push_back( qclab.qgates.U3(1, gammaL2, betaL2, deltaL2 ) ) ; 
        circuit0.push_back( qclab.qgates.CNOT(1,0) ) ; 
        circuit0.push_back( qclab.qgates.RotationZ(0, phiU ) ) ; 
        circuit0.push_back( qclab.qgates.RotationY(1, psiU ) ) ; 
        circuit0.push_back( qclab.qgates.CNOT(0,1) ) ; 
        circuit0.push_back( qclab.qgates.RotationY(1, thetaU ) ) ;
        circuit0.push_back( qclab.qgates.CNOT(1,0) ) ; 
        circuit0.push_back( qclab.qgates.U3(0, gammaR1, betaR1, deltaR1 ) ) ; 
        circuit0.push_back( qclab.qgates.U3(1, gammaR2, betaR2, deltaR2 ) ) ; 
        
        % circuit_global_phase = qclab.QCircuit( circuit.nbQubits ) ; 
        % circuit_global_phase.push_back( qclab.qgates.RotationZ(0, global_phase_angle ) ) ; 
        % circuit_global_phase.push_back( circuit ) ; 
        circuit.push_back( circuit0 ) ; 
        % circuit = circuit_global_phase ; 
        
        Delta = eye(4) ; 

        if logging 
            CNOT_count = 3 ; 
        end

        % M = circuit0.matrix ;
        % norm( U * Delta0' .* global_phase0  - Delta' *  M(1:4,1:4) .* global_phase ) 
        
    end
    if ~logging, CNOT_count = false ; end
end





%% test uniformly_controlled_rotation_z()


% function [circuit, CNOT_count] = uniformly_controlled_rotation_z( circuit, D, n, mode, logging )
% % uniformly controlled rotation R_Z 
% % Input: 
% %   circuit -- qclab circuit 
% %   D       -- diagonal matrix 
% %   n       -- circuit size 
% %   mode    -- 'L', 'M', 'R' (Three kinds of uniformly controlled rotation R_Z in `Beyond quantum Shannon decomposition') 
% %   mode L repsent the circuit -Rz-CNOT-...-CNOT-Rz- (Rz first uniformly controlled rotation, the last CNOT has been merged) 
% %   mode M repsent the circuit -Rz-CNOT-...-CNOT-Rz-CNOT- (Rz first uniformly controlled rotation) 
% %   mode R repsent the circuit -Rz-CNOT-...-CNOT-Rz- (CNOT first uniformly controlled rotation, the first CNOT has been merged) 
% % ---------------------------------------------------
% % Output: 
% %   circuit -- qclab circuit 
% % ---------------------------------------------------
% % test program:
%     % circuit = qclab.QCircuit(3) ; 
%     % circuit = uniformly_controlled_rotation_z( circuit, DA, n, 'M' ) ; 
%     % norm( circuit.matrix - [DA, zeros(size(DA)); zeros(size(DA)), DA'] ) 
%     % 
%     % circuit = qclab.QCircuit(3) ; 
%     % circuit = uniformly_controlled_rotation_z( circuit, DA, n, 'L' ) ; 
%     % circuit.push_back( qclab.qgates.CNOT(1,0) ) ; 
%     % norm( circuit.matrix - [DA, zeros(size(DA)); zeros(size(DA)), DA'] ) 
%     % 
%     % circuit = qclab.QCircuit(3) ; 
%     % circuit.push_back( qclab.qgates.CNOT(1,0) ) ; 
%     % circuit = uniformly_controlled_rotation_z( circuit, DA, n, 'R' ) ; 
%     % norm( circuit.matrix - [DA, zeros(size(DA)); zeros(size(DA)), DA'] ) 
% % --------------------------------------------------- 
% 
%     % Step 1. Compute rotation angles
%     phi = UniformlyRotationAngle( -angle(diag(D)) ).*2 ; 
% 
%     % Step 2. Construct quantum circuit
%     circuit0 = qclab.QCircuit( n , circuit.nbQubits - n ) ; 
%     if mode == 'L' || mode == 'M' 
%         for k = 1 : size(phi,1) 
%             circuit0.push_back( qclab.qgates.RotationZ(0, phi(k)) ) ; 
%             if mode ~= 'L' || k ~= size(phi,1)
%                 circuit0.push_back(qclab.qgates.CNOT(ctrl_pos(k, n-1), 0)) ; 
%             end
%         end
%     else % mode == 'R'
%         for k = size(phi,1):-1:1
%             if k ~= size(phi,1)
%                 circuit0.push_back(qclab.qgates.CNOT(ctrl_pos(k, n-1), 0)) ; 
%             end
%             circuit0.push_back(qclab.qgates.RotationZ(0, phi(k))) ; 
%         end
%     end
%     circuit.push_back( circuit0 ) ; 
% 
%     if logging
%         if mode == 'L' || mode == 'R'
%             CNOT_count = size(phi,1) - 1 ;
%         else % mode == 'M'
%             CNOT_count = size(phi,1) ; 
%         end
%     end
% 
% end



% function ctrl = ctrl_pos(i, n)
%     ctrl = n - log2(bitxor(grayCode(i-1),grayCode(i))) ;
%     if i == pow2(n)
%         ctrl = 1;
%     end
% end % end of ctrl_pos


% function [alpha, gamma, beta, delta] = u2_decomposition(U)
%     % U2_DECOMPOSITION_ROBUST Decomposes a U(2) matrix into Rz(psi) * Ry(theta) * Rz(phi) with global phase
%     % Input: U - 2x2 unitary matrix
%     % Outputs: alpha - global phase
%     %          beta, delta, gamma - rotation angles for Ry and Rz gates
%     % such that 
%     % ---------------------------------------------------------------------
%     % circuit = qclab.QCircuit(1)
%     % circuit.push_back( qclab.qgates.U3( 0, gamma, beta, delta ) )
%     % norm( exp(1j*alpha) .* circuit.matrix - U ) 
%     % ---------------------------------------------------------------------
%     % Step 1. Compute parameter 
%     % such that 
%     % U = exp(1i*alpha) * [exp(-1i*beta/2), 0; 0, exp(1i*beta/2)] * [cos(gamma/2), -sin(gamma/2); sin(gamma/2), cos(gamma/2)] * [exp(-1i*delta/2), 0; 0, exp(1i*delta/2)];
%     alpha = angle(det(U)) / 2;
%     beta = (angle( -U(2,1) / U(1,2) ) + angle( U(2,2) / U(1,1) )) ./ 2 ;
%     delta = (-angle( -U(2,1) / U(1,2) ) + angle( U(2,2) / U(1,1) )) ./ 2  ;
%     gamma = angle(exp(1i*(-alpha+beta/2+delta/2)).*U(1,1) + 1j.*exp(1i*(-alpha-beta/2+delta/2)).*U(2,1)).*2;
% 
%     % Step 2. QCLAB U3 function global phase shift
%     alpha = alpha - beta / 2 - delta / 2 ;
% 
% end







% function [thetat] = UniformlyRotationAngle(theta)
% % Compute the uniformly controlled rotation
% % thetat = (M^n)^(-1)*theta
%     thetat = grayPermutation( sfwht( theta ) ) ; 
% end % end of UniformlyRotationAngle
% 
% function [ b ] = grayPermutation( a )
%   k = log2( size(a, 1) ) ; 
%   b = zeros( size(a) );
%   for i = 0 : 2^k - 1
%     b( i + 1, : ) = a( grayCode( i ) + 1, : );
%   end
% end % end of grayPermutation
% 
% function [ a ] = sfwht( a )
% % Scaled fast Walsh-Hadamard transform
%   k = log2(size(a, 1) ) ;
%   for h = 1:k
%     for i = 1:2^h:2^k
%       for j = i:i+2^(h-1)-1
%         x = a( j ,: );
%         y = a( j + 2^(h-1) ,: );
%         a( j ,: ) = ( x + y ) / 2 ;
%         a( j + 2^(h-1) ,: ) = ( x - y ) / 2 ;
%       end
%     end
%   end
% end % end of sfwht
% 
% function x = grayCode(x) 
%     x = bitxor(x,bitshift(x,-1)) ; 
% end % end of grayCode 
