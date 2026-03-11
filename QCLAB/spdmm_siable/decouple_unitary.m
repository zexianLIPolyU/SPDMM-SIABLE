function [U,D,V] = decouple_unitary(U1,U2)
% Decoupling diagonal unitary
% (U1 \oplus U2) = (U \oplus U) * (D \oplus D^\dagger) * (V \oplus V) 

    [U, D] = eig( U1*U2' ) ;
    D = diag(sqrt(diag(D))) ;
    V = D * U' * U2 ;
    % V = U2' * U * D' ; 

end
