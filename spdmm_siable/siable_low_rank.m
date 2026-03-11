function [circuit, subnormalization, info] = siable_low_rank(Matrix, rank, logging, offset)
% SIABLE -- Single Ancilla Block Encodings for low-rank matrix
% INPUT
% -----
% Matrix:       matrix to be block encoded 
% logging:      true/false, if true info will log information about compression
% offset:       starting position in the circuit. The default value for 'offset' is 0.
%
% OUTPUT
% ------
% circuit:                  QCLAB circuit that block encodes A 
% subnormalization:         subnormalization factor of this block-encoding 
% CNOT_count:               CNOT's number 
%
%
% Reference: 
% Improving C-NOT Counts for Quantum State Preparation and Block Encoding via Diagonal Matrix Migration
% Minimal universal two-qubit controlled-NOT-based circuits 
% Beyond quantum Shannon decomposition: Circuit construction for n-qubit gates based on block-ZXZ decomposition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin <= 2 
        offset = 0 ; 
    end 
    if nargin <= 1 
        logging = false ; 
    end 
    if logging 
        info.nCNOT = 0 ; 
    else 
        info = false ; 
    end 


    % Step 1. subnormalization 
    subnormalization = norm(Matrix, 2) ; 
    Matrix = Matrix ./ subnormalization ; 
    N = size(Matrix, 1) ; 
    n = log2( N ) ; 
    assert( rank <= pow2(n-1) ) ; 
    assert( n >= 2 ) ; 
    
    assert( N == pow2(n) ) ; 
    assert( N == size(Matrix, 2) ) ; 
    
    % Matrix = Matrix ./ subnormalization ; 
    % Step 2 Compute the circuit parameters 
    % Find two unitaries such that norm((U1+U2)./2 - Matrix) 
    A1 = Matrix ; 
    [W, S, V] = svd(A1) ; 
    S1 = diag(S) ; 
    S1 = S1(1:rank) ; 
    if rank ~= pow2(n) 
        S1(pow2(n)) = 0 ; 
    end 
    D = exp( 1j.*acos(S1) ) ; 
    
    % U1 = W * diag(D) * V' ; 
    % U2 = W * diag(D') * V' ; 
    % % Decouple unitary
    % [U, D, V] = decouple_unitary(U1, U2) ; 
    % 
    % % Step 3 Construct quantum circuit 
    NumQubits = n + 1 ; 
    circuit = qclab.QCircuit( NumQubits, offset ) ; 
    
    % [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit( circuit, V, n, 1, eye(4), false, logging ) ;
    % if logging, CNOT_count = CNOT_count0 + CNOT_count ; end

    up_to_diagonal_matrix = true ; 
    [circuitL, global_phase0, info0, Delta] = ciruict_isometry( V, rank, n, logging, up_to_diagonal_matrix ) ; 
    if logging, info.nCNOT = info.nCNOT + info0.nCNOT ; end 
    % Delta - diag( circuitL.matrix * V )
    % M = circuitL.matrix;  diag( Delta' ) * M .*global_phase0 - V' 
    
    circuitM = qclab.QCircuit(n + 1, offset) ; 
    [circuitM, nCNOT0] = uniformly_controlled_rotation_z( circuitM, D, n+1, 'M', logging ) ; 
    if logging, info.nCNOT = info.nCNOT + nCNOT0 ; end 
    
    up_to_diagonal_matrix = false ; 
    [circuitR, global_phase1, info0, ~] = ciruict_isometry( W * diag(Delta'), rank, n, logging, up_to_diagonal_matrix ) ; 
    if logging, info.nCNOT = info.nCNOT + info0.nCNOT ; end 
    % M = circuitR.matrix; global_phase1'.* M'  - W * diag(Delta') 

    circuit.push_back( qclab.qgates.RotationZ(0, -angle( global_phase0 * global_phase1' ).*2 ) ) ; 
    circuit.push_back( qclab.qgates.Hadamard(0) ) ; 
    circuit.push_back( circuitL ) ; 
    circuit.push_back( circuitM ) ; 
    circuit.push_back( circuitR.ctranspose() ) ; 

    % [circuit, global_phase, ~, CNOT_count0] = recursive_decouple_quantum_circuit( circuit, U, n, global_phase, Delta, true, logging ) ; 
    % if logging, CNOT_count = CNOT_count0 + CNOT_count ; end
    circuit.push_back( qclab.qgates.Hadamard(0) ) ; 
    % circuit.push_back( qclab.qgates.RotationZ(0, -2.*angle(global_phase))
    % ) ; % global phase 

    % M =  circuit.matrix ; 
    % 
    % lowrank_approximation = W * diag(S1) * V' ; 
    % 
    % norm( lowrank_approximation - M(1:N, 1:N) ) 
    % 
    % % test programe 
    % M1 = kron([1,1;1,-1]./sqrt(2), W) * [diag(D),zeros(length(D));zeros(length(D)),diag(D')] * kron([1,1;1,-1]./sqrt(2),V') ;
    % norm( M1(1:N,1:N) - W * diag(S1) * V' )
    % 
    %  a = 0 ; 
    
end


function [circuit, global_phase, info, Delta] = ciruict_isometry( V, k, n, logging, up_to_diagonal_matrix )
% Generate a quantum circuit whose first k rows are the conjugate transport
% of the first k columns of V (up to a diagonal matrix) . 
% input: 
% V -- isometry size of pow2(n)*pow2(n) 
% k -- the column number to be generate 
% n -- the number of qubit 
% logging -- record the number of CNOT gates or not 
% up_to_diagonal_matrix -- generate the isometry up to a diagonal matrix or
% nopt.
% -------------------------------------------------------------------------
% output: 
% circuit -- the quantum circuit 
% global_phase -- the global phase 
% info -- the number of CNOT gates  
% -------------------------------------------------------------------------
% test program
%                                                                                                                                             
% M = circuitL.matrix ; 
% M * diag(Delta) .* global_phase0 - V' 
% the first few rows are zeros' row .
% -------------------------------------------------------------------------

    offset = 1 ; 
    circuit = qclab.QCircuit(n, offset) ; 
    if logging 
        info.nCNOT = 0 ; 
    else
        info = false ; 
    end 
    assert( k >= 1) ; 
    k = k - 1 ;
    for i = 0 : k  
        state0 = circuit.matrix * V( :, i + 1 ) ; 
        state = state0 ; 
        % global_phase = 1 ; 
        if i == 0 
            [circuit, global_phase, nCNOT0] = state_preparation( state0, offset, logging ) ; 
            if logging, info.nCNOT = info.nCNOT + nCNOT0 ; end 
            circuit = circuit.ctranspose() ; 
        else 
            Gk = qclab.QCircuit(n) ; 
            for s = 0 : n - 2 
                U2_sequence = decouple_sequence(i, s, state) ; 
                [ Gk, global_phase0, ~, nCNOT0 ] = reduced_multiplxor( n, U2_sequence, 0:n-2-s, n-s-1, Gk ) ; 
                if logging, info.nCNOT = info.nCNOT + nCNOT0 ; end 
                global_phase = global_phase * global_phase0 ;  
                state = Gk.matrix * state0 ; % state = global_phase * Gk.matrix * state0 ;
            end 
            corresponding_arrays = generate_state_index(i, n) ; % generate the non-zero index in a state 
            MCG_count = 0 ; 
            for s = 1 : size(corresponding_arrays, 2) 
                if mod( log2( bitxor(i, corresponding_arrays(s) ) ), 1 ) == 0 
                    [controls, target, controlState] = decouple_control_pair2(i, corresponding_arrays(s), n) ; 
                    Gk = MCG( target, decouple_onebit(0, [state(i + 1),state(corresponding_arrays(s) + 1)]), controls, controlState, Gk ) ; 
                    if logging
                        if n >= 6
                            % Reference: PHYS. REV. APPLIED 24, 044030 (2025)
                            info.nCNOT = info.nCNOT + 12*n - 32 ; 
                        else
                            % Reference: Circuit decomposition of multicontrolled special unitary single-qubit gates. Trans. Comp.-Aided Des. Integ. Cir. Sys., 43(3):802–811, March 2024.
                            if mod(n,2) % n is an odd number 
                                info.nCNOT = info.nCNOT + 20*n - 38 ; 
                            else % n is an even number
                                info.nCNOT = info.nCNOT + 20*n - 42 ; 
                            end 
                        end
                    end 
                    state = Gk.matrix * state0 ; % state = global_phase * Gk.matrix * state0 ; 

                    MCG_count = MCG_count + 1 ; 
                end 
            end 
            circuit.push_back(Gk) ; 
        end 
    end 
    

    Delta = diag( circuit.matrix * V ) ; 
    if ~up_to_diagonal_matrix
        phase = Delta( 1:k+1 ) ; 
        if mod(log2(k+1),1) ~= 0 
            phase(pow2(ceil(log2(k+1)))) = 0 ; 
        end 
        if k == 0 
            phase(2,1) = 0 ; 
        end 
        
        state_phase = angle(phase) ; 
        N = length( state_phase ) ; 
        m = log2( length(state_phase) ) ; 
        circuit_sim = true ; 
        HatTheta = AngleCompute('RZ', state_phase) ; % HatTheta for rotation-Z
        HatTheta(1:N) = HatTheta([2:N,1]) ;  
        % put the gloabl phase in the end  
        for i = 2 : m 
            row = pow2(i-1) : pow2(i)-1 ; 
            % Varphi(row) = UniformlyRotationAngle(Varphi(row));
            if ~isreal(state) 
                HatTheta(row) = UniformlyRotationAngle(HatTheta(row)) ; 
            end 
        end 
        
        % nRZ = 0 ; 
        circuit2 = qclab.QCircuit( m, n - m  ) ; 
        % circuit2.push_back(qclab.qgates.RotationZ(0,HatTheta(N))) ; 
        circuit2.push_back(qclab.qgates.RotationZ( 0, HatTheta(1) )) ; 
        % if logging, nRZ = nRZ + 1 ; end 
        for k = 1 : m-1 
            para_row_index = pow2(k):pow2(k+1)-1 ; 
            [circuit2, info0 ] = UniformRotation2( circuit2, 'RZ', HatTheta(para_row_index)', setdiff(0:m-1, k), k, logging, circuit_sim) ; 
            if logging 
                info.nCNOT = info.nCNOT + info0.nCNOT ; 
            end 
        end 
        circuit.push_back(circuit2.ctranspose) ; 
        global_phase = exp(1i./2*HatTheta(N)) ; 
        Delta = ones(pow2(n), 1) ; 
    else
        Delta(k+2:end) = ones(pow2(n)-k-1,1) ;
        global_phase = 1 ; 
    end 
end 




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
    
    % step 0. Compute the special case of [1,1j,1j,...,1j];

    % special_vec = ones(size(D)) .* 1j ; special_vec(1) = 1 ; 
    % if norm( D - special_vec) <= 1e-4 
    %     % Step 1. Compute rotation angles 
    %     phi = UniformlyRotationAngle( -angle(D) ).*2 ; 
    %     circuit0 = qclab.QCircuit( n , circuit.nbQubits - n ) ; 
    %     for k = 1 : length(phi) 
    %         circuit0.push_back( qclab.qgates.RotationZ(0, phi(length(phi) - k + 1)) ) ; 
    %         if mode ~= 'L' || k ~= length(phi) 
    %             circuit0.push_back( qclab.qgates.CNOT(ctrl_pos(k, n-1), 0) ) ; 
    %         end 
    %     end 
    %     circuit.push_back( circuit0 ) ; 
    % else


    % Step 1. Compute rotation angles 
    if size(D,1) < size(D,2) 
        D = D.' ; 
    end
    phi = UniformlyRotationAngle( -angle(D) ).*2 ; 

    % [ circuit, info, ~ ] = UniformRotation2( circuit, 'RZ', phi, 1:n-1, 0, logging, true ) ; 
    % CNOT_count = info.nCNOT ;
    % Step 2. Construct quantum circuit 
    circuit0 = qclab.QCircuit( n , circuit.nbQubits - n ) ; 
    if mode == 'L' || mode == 'M' 
        for k = 1 : length(phi) 
            circuit0.push_back( qclab.qgates.RotationZ(0, phi(k)) ) ; 
            if mode ~= 'L' || k ~= length(phi) 
                circuit0.push_back( qclab.qgates.CNOT(ctrl_pos(k, n-1), 0) ) ; 
            end 
        end 
    else % mode == 'R' 
        for k = length(phi):-1:1 
            if k ~= length(phi) 
                circuit0.push_back(qclab.qgates.CNOT(ctrl_pos(k, n-1), 0)) ; 
            end 
            circuit0.push_back(qclab.qgates.RotationZ(0, phi(k))) ; 
        end 
    end 
    circuit.push_back( circuit0 ) ; 


    if logging 
        if mode == 'L' || mode == 'R' 
            CNOT_count = length(phi) - 1 ; 
        else % mode == 'M'
            CNOT_count = length(phi) ; 
        end 
    end 

end 

%% UniformRotation2()
function [ circuit, info, parity_check ] = UniformRotation2( circuit, ctrl_type, para_seq, ctrl_index, targ_index, logging, circuit_sim ) 
% Input:    circuit     --  generated by qclab.QCircuit; 
%           ctrl_type   --  'RY' for Rotation-Y/ 'RZ' for Rotation-Z 
%           para_seq    --  parameter generated by Walsh-Hadamard transform 
%           ctrl_index  --  a vector contain the control 
%           logging     --  true/false, if true info will log information about compression 
% Output:   circuit     --  QCLAB circuit that block encodes A    
%           info        --  struct containing some info on compression algorithm and circuit

    n_count = log(size(para_seq, 2)) / log(2) ; 
    if strcmp( ctrl_type, 'RY' ) 
        G = @qclab.qgates.RotationY ; 
    elseif strcmp( ctrl_type, 'RZ' ) 
        G = @qclab.qgates.RotationZ ; 
    end 
    
    nG = 0 ; nCNOT = 0 ; 
    i = 1 ; 
    parity_check = int32(0) ; 
    while i <= pow2(n_count) 
        if any(para_seq(i) ~= 0) 
            % Add CNOTs based on parity_check
            [ circuit, num_CNOT ] = make_CNOT( circuit, parity_check, ctrl_index, targ_index, circuit_sim ) ; 
            nCNOT = nCNOT + num_CNOT ; 

            if circuit_sim 
                circuit.push_back( G(targ_index, para_seq(i)) ) ; 
            end 
            nG = nG + 1 ; 
            ctrl = ctrl_pos( i, n_count ) ; 
            % update parity check
            parity_check = bitset( int32(0), ctrl_index(ctrl) + 1, int32(1) ) ; 
            i = i + 1 ; 
        else 
            % update parity check
            while i <= pow2(n_count) && all( para_seq(:,i) == 0 ) 
                ctrl = ctrl_pos(i, n_count) ; 
                if bitget( parity_check, ctrl_index(ctrl) + 1 ) 
                    parity_check = bitset( parity_check, ctrl_index(ctrl) + 1, int32(0) ) ; 
                else 
                    parity_check = bitset( parity_check, ctrl_index(ctrl) + 1, int32(1) ) ; 
                end 
                i = i + 1 ; 
            end 
        end 
    end 
    % Add CNOTs based on parity_check in the final
    [ circuit, num_CNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim ) ; 
    nCNOT = nCNOT + num_CNOT ; 
    if logging 
        info = struct() ; 
        info.nG = nG ; 
        info.nCNOT = nCNOT ; 
    else 
        info = false ; 
    end 
end % end of UniformRotation2



function ctrl = ctrl_pos(i, n) 
    ctrl = n - log2(bitxor(grayCode(i-1),grayCode(i))) ; 
    if i == pow2(n) 
        ctrl = 1 ; 
    end 
end % end of ctrl_pos



function [correspoinding_arrays] = generate_state_index(i, n) 
    correspoinding_arrays = zeros(1, max(1,ceil(log2(pow2(n-1)-i))) ) ; 
    k = 1 ; 
    for t = 1 : size(correspoinding_arrays, 2) 
        if t == 1 
            correspoinding_arrays(t) = i + pow2(n-t) ; 
        else 
            while i + pow2(n-k) >= pow2(n-1) 
                k = k + 1 ; 
            end 
            correspoinding_arrays(t) = i + pow2(n-k) ; 
            k = k + 1 ; 
        end 
    end 
end 



function U2_sequence = decouple_sequence(i, s, state) 
% generate U(2) sequence whose uniformly controlled rotations decouple the
% s-th bit of state 

    n = log2(size(state,1)) ; 
    binary_char = dec2bin(i, n) ; 
    sequence_index = 0:pow2(s):pow2(n) - 1 ; 
    if s >= 1 
        sequence_index = sequence_index + ... 
            bin2dec(binary_char(end - min(s,size(binary_char,2)) + 1:end)) ; 
    end 
    sequence_index = sequence_index(sequence_index <= pow2(n) - 1) ; 
    U2_sequence = [] ; 
    if ceil(i/pow2(s+1)) >= 1 
        for t = 1 : ceil(i/pow2(s+1)) 
            U2_sequence = cat(3, U2_sequence, eye(2)) ; 
        end 
    end 
    for t = ceil(i/pow2(s+1)) + 1 : pow2(n - s - 1) 
        U2_sequence = cat(3, U2_sequence, decouple_onebit(mode(bin2dec(binary_char(n-s)),2), [state(sequence_index(2*t-1) + 1),state(sequence_index(2*t) + 1)])) ; 
    end 
end 



function [controls, target, controlState] = decouple_control_pair2(index1, index2, n) 
% Generate the controls index and their controlState 
    s = log2( bitxor(index1, index2) ) ; 
    assert( mod(s,1) == 0 ) ; 
    assert( s <= (n - 1) && s >= 1 ) ; 
    target = n - s - 1 ; 
    if s == n - 1 
        controls = 1 : n - 1 ; 
    else 
        controls = [ 0 : n - s - 2, n - s : n - 1 ] ; 
    end 
    binary_char = dec2bin(index1, n) ; 
    controlState = zeros(1 , n - 1) ; 
    for t = 1 : s 
        controlState( n - t ) = bin2dec(binary_char(n - t + 1)) ; 
    end 
    if s + 2 <= n 
        for t = s + 2 : n 
            controlState( n - t + 1 ) = bin2dec(binary_char(n - t + 1)) ; 
        end 
    end 
end 


function u2 = decouple_onebit(is_odd, onequbit_state) 
% Unitary to decouple one-qubit state (idea is derived from Lemma 2 in
% Quantum Circuits for Isometries) 
    if is_odd % odd 
        u2 = [onequbit_state(2), -onequbit_state(1); onequbit_state(1)', onequbit_state(2)'] ; 
    else % even 
        u2 = [onequbit_state(1)', onequbit_state(2)'; -onequbit_state(2), onequbit_state(1)] ; 
    end 
    u2 = u2 ./ norm(u2) ; 
end

function [ circuit, nCNOT ] = make_CNOT(circuit, parity_check, ctrl_index, targ_index, circuit_sim )  
% Add CNOTs based on parity_check in the final 
    if nargin <= 3 
        targ_index = 0 ; 
    end
    nCNOT = 0 ; 
    if parity_check ~= int32(0) 
        for j = 1 : numel(ctrl_index) 
            if bitget( parity_check, ctrl_index(j) + 1 ) 
                if circuit_sim 
                    circuit.push_back( qclab.qgates.CNOT( ctrl_index(j), targ_index ) ) ; 
                end
                nCNOT = nCNOT + 1 ; 
            end 
        end 
    end 
end % end of make_CNOT 








function [thetat] = UniformlyRotationAngle(theta)
% Compute the uniformly controlled rotation
% thetat = (M^n)^(-1)*theta
    thetat = grayPermutation( sfwht( theta ) ) ;
end % end of UniformlyRotationAngle

function [ b ] = grayPermutation( a )
  k = log2( size(a, 1) ) ; 
  b = zeros( 2^k, 1 );
  for i = 0 : 2^k - 1
    b( i + 1 ) = a( grayCode( i ) + 1 );
  end
end % end of grayPermutation

function [ a ] = sfwht( a )
% Scaled fast Walsh-Hadamard transform
  k = log2(size(a, 1) ) ;
  for h = 1:k
    for i = 1:2^h:2^k
      for j = i:i+2^(h-1)-1
        x = a( j );
        y = a( j + 2^(h-1) );
        a( j ) = ( x + y ) / 2 ;
        a( j + 2^(h-1) ) = ( x - y ) / 2 ;
      end
    end
  end
end % end of sfwht


function x = grayCode(x)
    x = bitxor(x,bitshift(x,-1));
end


