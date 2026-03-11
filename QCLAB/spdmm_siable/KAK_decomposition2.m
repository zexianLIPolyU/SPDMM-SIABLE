function [ L1, L2, R1, R2, theta, phi, psi, global_phase ] = KAK_decomposition2(U)
% Decouple U4 gate
% Using 2 CNOT gates to decompose a unitary in U(2) and a global phase
% Reference: Minimal universal two-qubit controlled-NOT-based circuits,
% Proposition 4 and Theorem 3.
% --------------------------------------
%   Input:
%       U -- unitary in U(4)
% --------------------------------------
%   Output:
%       L1,L2,R1,R2 -- SU(2);
%       Delta = CNOT12 * kron(eye(2),Rz(psi)) * CNOT12 ;
%       theta = ( eigs_U1(1) + eigs_U1(2) ) ./ 2 ;
%       phi = ( eigs_U1(1) - eigs_U1(2) ) ./ 2 ;
%   such that 
%
%   diag([exp(-1i.*psi./2),exp(1i.*psi./2),exp(1i.*psi./2),exp(-1i.*psi./2)] * U .* global_phase 
%   = kron(L1,L2) * CNOT12 * kron( Rx(theta), Rz(phi) ) * CNOT12 * kron(R1,R2)  
%
%
%   where Delta =
%   CNOT12 * kron(eye(2),Rz(psi)) * CNOT12 =
%   =  diag([exp(-1i.*psi./2),exp(1i.*psi./2),exp(1i.*psi./2),exp(-1i.*psi./2)]) .
%   Besides,
%       t = diag(gamma(U1.').');
%       psi = atan(imag(t(1) + t(2) + t(3) + t(4)) ./ real( t(1) - t(2) - t(3) + t(4)) ) ;
%   where 
%       norm( (t(1) + t(4)) * exp(-1i*psi) + (t(2) + t(3)) * exp(1i*psi) - trace(gamma(U1*Delta)) ),
%       eigs_U1 = sort(angle(eig( gamma(U1*Delta) )), 'descend') ;
% --------------------------------------
% test program:
    % [ L1, L2, R1, R2, theta, phi, psi, global_phase ] = KAK_decomposition4(U) ;
    % norm( kron( R1,R2) * CNOT12 * kron( Rx(theta), Rz(phi) ) * CNOT12 * kron(L1,L2) - diag([exp(-1i.*psi./2),exp(1i.*psi./2),exp(1i.*psi./2),exp(-1i.*psi./2)]) * U .* global_phase  )
% --------------------------------------


    if norm(det(U) - 1) < 1e-5
        global_phase = 1 ; 
    else
        global_phase = 1./det(U)^(1/4) ; 
    end
    U2 = U .* global_phase ; 
    % U1 = U2' ; 

    % Pauli matrices
    X = [0 1; 1 0] ; 
    Y = [0 -1i; 1i 0] ; 
    CNOT12 = kron([0 0; 0 1], X) + kron([1 0; 0 0], eye(2)) ; 
    Rz = @(x) [exp(-1i*x./2), 0; 0, exp(1i*x./2)] ; 
    Rx = @(x) expm(-1i*x/2.*X) ; 
    E = [1 1j 0 0; 0 0 1j 1; 0 0 1j -1; 1 -1j 0 0] ./ sqrt(2) ; 
    gamma = @(u) u * kron(Y, Y) * u.' * kron(Y,Y) ; 

    % Algorithm 
    t = diag(gamma(U2)) ; 
    psi = atan(imag(t(1) + t(2) + t(3) + t(4)) ./ real( t(1) - t(2) - t(3) + t(4)) ) ; 
    Delta = CNOT12 * kron(eye(2),Rz(psi)) * CNOT12 ; 
    % norm( (t(1) + t(4)) * exp(-1i*psi) + (t(2) + t(3)) * exp(1i*psi) - trace(gamma(U1*Delta)) )
    eigs_DeltaU = sort(angle(eig( gamma(Delta*U2) )), 'descend') ; 
    theta = ( eigs_DeltaU(1) + eigs_DeltaU(2) ) ./ 2 ; 
    phi = ( eigs_DeltaU(1) - eigs_DeltaU(2) ) ./ 2 ; 
    w = CNOT12 * kron( Rx(theta), Rz(phi) ) * CNOT12 ; 
    u = E \ Delta * U2 * E; 
    v = E \ w * E; 
    [ua,va] = eig(u*u.') ; 
    [ub,vb] = eig(v*v.') ; 
    [~, order_a] = sort(angle(diag(va)), 'descend') ; 
    [~, order_b] = sort(angle(diag(vb)), 'descend') ; 
    ua = real(ua(:,order_a)) ; 
    ub = real(ub(:,order_b)) ; 
    if norm(det(ua) + 1) < 1e-5, ua(:,1) = ua(:,1) .* -1; end 
    if norm(det(ub) + 1) < 1e-5, ub(:,1) = ub(:,1) .* -1; end 
    c = real(v'*ub.'*ua.'*u) ; 
    if norm(det(c) + 1) < 1e-5, c(:,1) = c(:,1) .* -1; end 
    R = E*ua*ub/E ; 
    L = E*c/E ; 
    [L1, L2] = decompose_local_product(L) ;
    [R1, R2] = decompose_local_product(R) ;
    % test programm:
    % norm( kron( R1,R2) * CNOT12 * kron( Rx(theta), Rz(phi) ) * CNOT12 * kron(L1,L2) - Delta * U2  )


end

function [A, B] = decompose_local_product(U) 
% decompose two-qubit product $U\in SU(2)\otimes SU(2)$
% such that $U = kron( A, B )$, where $A, B \in SU(2)$.
    if det(U(1:2,1:2)) >= 1e-4
        a_11 = sqrt(det(U(1:2,1:2))) ; 
        B = U(1:2,1:2) ./ a_11 ; 
        a_12 = trace(U(1:2,3:4)*B')/2 ; 
        a_21 = trace(U(3:4,1:2)*B')/2 ; 
        a_22 = trace(U(3:4,3:4)*B')/2 ; 
        A = [a_11, a_12; a_21, a_22] ;
    else
        a_12 = sqrt(det(U(1:2,3:4))) ; 
        B = U(1:2,3:4) ./ a_12 ; 
        a_11 = trace(U(1:2,1:2)*B')/2 ; 
        a_21 = trace(U(3:4,1:2)*B')/2 ; 
        a_22 = trace(U(3:4,3:4)*B')/2 ; 
        A = [a_11, a_12; a_21, a_22] ; 
    end

end