function [circuit, CNOT_count] = uniformly_controlled_rotation_z( circuit, D, n, mode, logging )
% uniformly controlled rotation R_Z 
% Input: 
%   circuit -- qclab circuit 
%   D       -- diagonal matrix 
%   n       -- circuit size 
%   mode    -- 'L', 'M', 'R' (Three kinds of uniformly controlled rotation R_Z in `Beyond quantum Shannon decomposition') 
%   mode L repsent the circuit -Rz-CNOT-...-CNOT-Rz- (Rz first uniformly controlled rotation, the last CNOT has been merged) 
%   mode M repsent the circuit -Rz-CNOT-...-CNOT-Rz-CNOT- (Rz first uniformly controlled rotation) 
%   mode R repsent the circuit -Rz-CNOT-...-CNOT-Rz- (CNOT first uniformly controlled rotation, the first CNOT has been merged) 
% ---------------------------------------------------
% Output: 
%   circuit -- qclab circuit 
% ---------------------------------------------------
% test program:
    % circuit = qclab.QCircuit(3) ; 
    % circuit = uniformly_controlled_rotation_z( circuit, DA, n, 'M' ) ; 
    % norm( circuit.matrix - [DA, zeros(size(DA)); zeros(size(DA)), DA'] ) 
    % 
    % circuit = qclab.QCircuit(3) ; 
    % circuit = uniformly_controlled_rotation_z( circuit, DA, n, 'L' ) ; 
    % circuit.push_back( qclab.qgates.CNOT(1,0) ) ; 
    % norm( circuit.matrix - [DA, zeros(size(DA)); zeros(size(DA)), DA'] ) 
    % 
    % circuit = qclab.QCircuit(3) ; 
    % circuit.push_back( qclab.qgates.CNOT(1,0) ) ; 
    % circuit = uniformly_controlled_rotation_z( circuit, DA, n, 'R' ) ; 
    % norm( circuit.matrix - [DA, zeros(size(DA)); zeros(size(DA)), DA'] ) 
% --------------------------------------------------- 

    % Step 1. Compute rotation angles
    phi = UniformlyRotationAngle( -angle(diag(D)) ).*2 ; 

    % Step 2. Construct quantum circuit
    circuit0 = qclab.QCircuit( n , circuit.nbQubits - n ) ; 
    if mode == 'L' || mode == 'M' 
        for k = 1 : size(phi,1) 
            circuit0.push_back( qclab.qgates.RotationZ(0, phi(k)) ) ; 
            if mode ~= 'L' || k ~= size(phi,1)
                circuit0.push_back(qclab.qgates.CNOT(ctrl_pos(k, n-1), 0)) ; 
            end
        end
    else % mode == 'R'
        for k = size(phi,1):-1:1
            if k ~= size(phi,1)
                circuit0.push_back(qclab.qgates.CNOT(ctrl_pos(k, n-1), 0)) ; 
            end
            circuit0.push_back(qclab.qgates.RotationZ(0, phi(k))) ; 
        end
    end
    circuit.push_back( circuit0 ) ; 

    if logging
        if mode == 'L' || mode == 'R'
            CNOT_count = size(phi,1) - 1 ;
        else % mode == 'M'
            CNOT_count = size(phi,1) ; 
        end
    else
        CNOT_count = false ; 
    end

end



function ctrl = ctrl_pos(i, n)
    ctrl = n - log2(bitxor(grayCode(i-1),grayCode(i))) ;
    if i == pow2(n)
        ctrl = 1;
    end
end % end of ctrl_pos


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







function [thetat] = UniformlyRotationAngle(theta)
% Compute the uniformly controlled rotation
% thetat = (M^n)^(-1)*theta
    thetat = grayPermutation( sfwht( theta ) ) ; 
end % end of UniformlyRotationAngle

function [ b ] = grayPermutation( a )
  k = log2( size(a, 1) ) ; 
  b = zeros( size(a) );
  for i = 0 : 2^k - 1
    b( i + 1, : ) = a( grayCode( i ) + 1, : );
  end
end % end of grayPermutation

function [ a ] = sfwht( a )
% Scaled fast Walsh-Hadamard transform
  k = log2(size(a, 1) ) ;
  for h = 1:k
    for i = 1:2^h:2^k
      for j = i:i+2^(h-1)-1
        x = a( j ,: );
        y = a( j + 2^(h-1) ,: );
        a( j ,: ) = ( x + y ) / 2 ;
        a( j + 2^(h-1) ,: ) = ( x - y ) / 2 ;
      end
    end
  end
end % end of sfwht

function x = grayCode(x) 
    x = bitxor(x,bitshift(x,-1)) ; 
end % end of grayCode 
