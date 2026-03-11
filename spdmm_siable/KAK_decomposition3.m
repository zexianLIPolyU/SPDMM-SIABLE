function [ L1, L2, R1, R2, phi, psi, theta, global_phase ] = KAK_decomposition3(U0)
% ----------------------------------------------------------
% test programme:
    % norm( global_phase .* U0 - kron(R1, R2) * CNOT21 * kron(eye(2), Ry(theta)) * CNOT12 * kron( Rz(phi), Ry(psi) ) * CNOT21 * kron(L1,L2) )
% ----------------------------------------------------------

    % Step 0. Map U(4) to SU(4) (and phase)
    SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
    U = SWAP * U0 ; 
    global_phase = 1./ det(U)^(1/4) ; 
    U = U .* global_phase ; 

    % Step 1. Unconjugate U into the magic basis M
    M = [1 0 0 1j; 0 1j 1 0; 0 1j -1 0; 1 0 0 -1j] ./ sqrt(2);
    U_prime = M'*U*M ;
    
    % Step 2. Isolating the maximal torus
    M_squared = U_prime.'*U_prime ;

    % Step 3. Extracting $K_2$
    [P, D] = eig(M_squared) ;

    if norm(det(P) + 1) <= 1e-5
        P(:,1) = -1.*P(:,1) ;
    end
    K2 = P' ;

    % Step 4: Extracting $A$
    A = sqrt(diag(D)) ;

    if norm(prod(A) + 1) <= 1e-5
        A(1,1) = A(1,1) .* -1;
    end

    % Step 5: Extracting $K_1$
    K1 = U_prime * K2' * diag(A)' ;

    % Step 6: Extracting Local Gates
    R = M*K1*M' ;
    L = M*K2*M' ;
    [R2, R1] = decompose_local_product(R) ;
    [L1, L2] = decompose_local_product(L) ;

    % Step 7: Extracting the Canoncial Parameters
    theta_vec = angle(A) ;
    c1 = theta_vec(1) + theta_vec(2) ;
    c2 = -theta_vec(1) - theta_vec(3) ;
    c3 = -theta_vec(2) - theta_vec(3) ;
    % norm(U - global_phase.* kron(L1,L2) * expm(1/2.*1j.*(c1.*kron(X,X)+c2.*kron(Y,Y)+c3.*kron(Z,Z))) * kron(R1,R2) )
    
    % Step 8: Extract rotation parameters 
    X = [0 1; 1 0] ; 
    Y = [0 -1i; 1i 0] ; 
    Z = [1 0; 0 -1] ; 
    E = [1 1j 0 0; 0 0 1j 1; 0 0 1j -1; 1 -1j 0 0] ./ sqrt(2) ; 
    Sz = expm(-1j.*[1,0;0,-1].*pi./4) ; 
    Sz_s = expm(1j.*[1,0;0,-1].*pi./4) ; 
    Ud = (E' * expm(1/2.*1j.*(c1.*kron(X,X)+c2.*kron(Y,Y)+c3.*kron(Z,Z))) * E) ; % diagonal matrix
    angles = angle(diag(Ud)) ;
    
    R1 = R1 * Sz_s ; 
    L2 = Sz * L2 ; 
    phi = - angles(1) - angles(2) ; 
    psi = - angles(1) - angles(3) ; 
    theta = angles(2) + angles(3) ; 
    
end

function [A, B] = decompose_local_product(U) 
% decompose two-qubit product $U\in SU(2)\otimes SU(2)$
% such that $U =kron( A, B )$, where $A, B \in SU(2)$.
    a_11 = sqrt(det(U(1:2,1:2))) ;
    B = U(1:2,1:2) ./ a_11 ;
    a_12 = trace(U(1:2,3:4)*B')/2 ;
    a_21 = trace(U(3:4,1:2)*B')/2 ;
    a_22 = trace(U(3:4,3:4)*B')/2 ;
    A = [a_11, a_12; a_21, a_22] ;
end