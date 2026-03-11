function [circuit, global_phase, CNOT_count] = state_preparation( vector, offset, logging )
% state_preparation 
%   reference:  
%   0. Improving C-NOT Counts for Quantum State Preparation and Block Encoding via Diagonal Matrix Migration
%   1. Quantum-state preparation with universal gate decompositions, Phys. Rev. A 83, 032302 (2011). 
%   2. Quantum Circuits for Isometries, Phys. Rev. A 93, 032318 (2016). 
% --------------------------------------------------
% test program:
% state_complex = randn(N,1) + randn(N,1) .* 1j ; 
% state_complex = state_complex ./ norm(state_complex,2) ;
% [circuit, global_phase, CNOT_count] = state_preparation( state_complex, logging ) ; 
% M = circuit.matrix ;
% circuit_state = M(:,1) .* global_phase ; 
% norm( circuit_state - state_complex )
% -------------------------------------------------- 
    n = log2( size( vector, 1 ) ) ; 
    if nargin <= 1 
        logging = false ; offset = 0 ; 
    elseif nargin <= 2 
        logging = false ; 
    end 
    if logging 
        CNOT_count = 0 ; 
    else 
        CNOT_count = false ; 
    end 
    global_phase = 1 ; 
    if n == 1 
        SU2 = [vector, [-vector(2)'; vector(1)']] ; 
        [ alphaSU2, gammaSU2, betaSU2, deltaSU2 ] = u2_decomposition( SU2 ) ;
        % [circuit, info] = binary_tree_statepreparation( vector, offset, logging, true ) ; 
        circuit = qclab.QCircuit( 1, offset ) ; 
        circuit.push_back( qclab.qgates.U3( 0, gammaSU2, betaSU2, deltaSU2 ) ) ; 
        global_phase = exp(1j.*alphaSU2) ; 
        % if logging, CNOT_count = CNOT_count + info.circ.nCNOT ; end 
    elseif n >= 2 
        reshape_vec = reshape( vector, pow2(floor(n/2)), [] ) ; 
        [U, S, V] = svd( reshape_vec.' ) ; 
        % norm( U*S*V' - reshape_vec.' )  

        NumQubits = n ; 
        circuit = qclab.QCircuit( NumQubits, offset ) ; 
        circuitUV = qclab.QCircuit( NumQubits, 0 ) ; 

        % Phase 2 -- COPY 
        for i = 1 : floor(n/2)  
            circuitUV.push_back( qclab.qgates.CNOT(i + ceil(n/2) - floor(n/2) - 1, i + ceil(n/2) - 1 ) ) ; 
        end 
        if logging, CNOT_count = CNOT_count + floor(n/2) ; end 

        % Phase 3 -- Prepare U 
        circuit0 = qclab.QCircuit( ceil(n/2), 0 ) ; 
        if ceil(n/2) == 1 
            [alpha, gamma, beta, delta ] = u2_decomposition(U') ; 
            circuit0.push_back( qclab.qgates.U3(0, gamma, beta, delta) ) ; 
            global_phase0 = exp( 1j.*alpha ) ; 
        else 
            if ceil(n/2) > floor(n/2) && ceil(n/2) > 2  % Prepare an unitary whose first pow2(floor(n/2)) columns are the isometry form pow2(ceil(n/2)) to pow2(floor(n/2)) 
                [circuit0, global_phase0, DeltaU, CNOT_count0] = isometry(U', circuit0, ceil(n/2), logging) ; 
                if logging, CNOT_count = CNOT_count + CNOT_count0 ; end 
                % C = circuit.ctranspose ; 
                % M = C.matrix * kron(eye(pow2(ceil(n/2)-2)), Delta) .* global_phase' ; 
                % norm( M(:,1:pow2(floor(n/2))) - U(:,1:pow2(floor(n/2))) ) 
            else % prepare a unitary form pow2(n/2) to pow2(n/2) 
                is_last_U = false ; 
                [circuit0, global_phase0, DeltaU, CNOT_count0] = recursive_decouple_quantum_circuit(circuit0, U', ceil(n/2), 1, eye(4), is_last_U, true ) ; 
                if logging, CNOT_count = CNOT_count + CNOT_count0 ; end 
                % C = circuit0.ctranspose ; norm(C.matrix * kron(eye(pow2(ceil(n/2)-2)), DeltaU) .* global_phase0' - U) 
            end
        end 
        global_phase = global_phase * global_phase0' ; 
        circuitUV.push_back( circuit0.ctranspose() ) ; 

        % Phase 4 -- Prepare conj(V) 
        circuit0 = qclab.QCircuit( floor(n/2), ceil(n/2) ) ; 
        if floor(n/2) == 1 
            [alpha, gamma, beta, delta ] = u2_decomposition( conj(V') ) ; 
            circuit0.push_back( qclab.qgates.U3(0, gamma, beta, delta) ) ; 
            global_phase0 = exp( 1j.*alpha ) ; 
        else 
            is_last_U = false ; 
            [circuit0, global_phase0, DeltaV, CNOT_count0] = recursive_decouple_quantum_circuit(circuit0, conj(V'), floor(n/2), 1, eye(4), is_last_U, true ) ; 
            if logging, CNOT_count = CNOT_count + CNOT_count0 ; end 
        end 
        global_phase = global_phase * global_phase0' ; 
        circuitUV.push_back( circuit0.ctranspose() ) ; 

        % Phase 1 -- recursive prepare 
        %   diag( kron(eye(pow2(n-2)), DeltaU') * S * kron(eye(pow2(n-2)), DeltaV) ) 
        if n == 2 
            vector = diag( S ) ; 
        elseif n == 3 
            vector = diag( kron(eye(pow2(ceil(n/2)-2)), DeltaU) * S ) ; 
        else % n >= 4 
            vector = diag( kron(eye(pow2(ceil(n/2)-2)), DeltaU) * S * kron(eye(pow2(floor(n/2)-2)), DeltaV) ) ; 
        end 
        [circuit0, global_phase0, CNOT_count0] = state_preparation( vector, ceil(n/2)-floor(n/2), logging ) ; 
        global_phase = global_phase * global_phase0 ; 
        if logging, CNOT_count = CNOT_count + CNOT_count0 ; end 
        circuit.push_back( circuit0 ) ; 

        % Phase 2-4 -- Prepare U, conj(V) 
        circuit.push_back( circuitUV ) ; 

    end

end



function [circuit, global_phase, Delta, CNOT_count] = isometry(U, circuit, n, logging) 
% Prepare an pow2(floor(n/2)) * pow2(floor(n/2)) unitary whose first pow2(floor(n/2)) rows are the isometry form pow2(n) to pow2(floor(n/2)) 
%   but attention that above n is different from the input of isometry() 
% Input: 
%   U            -- unitary to be decoupled, whose whose first pow2(floor(n/2)) columns are the isometry form pow2(ceil(n/2)) to pow2(floor(n/2)) 
%   n            -- parameter 
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict 
%   is_last_U    -- boolean logical that whether U is the last unitary in the decomposition 
% --------------------------------------------------- 
% Output: 
%   circuit      -- qclab circuit 
%   global_phase -- the global phase of the circuit 
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict 
%   CNOT_count   -- the CNOT count 
% ---------------------------------------------------
% test program: 
    % n = 5 ; 
    % U = generate_U2n( ceil(n/2) ) ; 
    % circuit = qclab.QCircuit( ceil(n/2) ); 
    % [circuit, global_phase, Delta, CNOT_count] = isometry( U', circuit, ceil(n/2), true ) ;  
    % C = circuit.ctranspose ; 
    % M = C.matrix * kron(eye(pow2(ceil(n/2)-2)), Delta) .* global_phase' ; 
    % norm( M(:,1:pow2(floor(n/2))) - U(:,1:pow2(floor(n/2))) ) 
% --------------------------------------------------- 

    assert( logging == true ) ; 
    if logging, CNOT_count = 0 ; end
    
    % Step 1. Compute the rotation parameters 
    [~, B1, B2, C] = transfomed_ZXZ(U) ; 
    % [A, B1, B2, C] = transfomed_ZXZ(U) ; 
    % U = [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), A]*kron(H,eye(pow2(n-1))) * [B1,zeros(pow2(n-1));zeros(pow2(n-1)),B2] * kron(H, eye(pow2(n-1))) * [ eye(pow2(n-1)), zeros(pow2(n-1));zeros(pow2(n-1)), C ] ;
    [WC, DC, VC] = decouple_unitary( eye(pow2(n-1)), C ) ; 
    % norm([WC, zeros(pow2(n-1)); zeros(pow2(n-1)), WC] * [DC, zeros(pow2(n-1)); zeros(pow2(n-1)), DC'] * [VC, zeros(pow2(n-1)); zeros(pow2(n-1)), VC] - [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), C] )  
    Z = [1 0; 0 -1] ; 
    [WB, DB, VB] = decouple_unitary( B1*WC, B2*WC*kron(Z,eye(pow2(n-2))) ) ; 
    % norm([WB, zeros(pow2(n-1)); zeros(pow2(n-1)), WL] * [DB, zeros(pow2(n-1)); zeros(pow2(n-1)), DB'] * [VB, zeros(pow2(n-1)); zeros(pow2(n-1)), VB] - [B1*WC, zeros(pow2(n-1)); zeros(pow2(n-1)), B2*WC*kron(Z,eye(pow2(n-2)))] )  
    
    % Step 2. Construct quantum circuit 
    global_phase = 1 ; Delta = eye(4) ; 
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
    [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, WB, n-1, global_phase, Delta, false, logging) ; 
    if logging, CNOT_count = CNOT_count0 + CNOT_count ; end
    circuit.push_back( qclab.qgates.Hadamard(circuit.nbQubits - n) ) ;

end 

