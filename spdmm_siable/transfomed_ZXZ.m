function [A, B1, B2, C] = transfomed_ZXZ(U)
% Transfomed block-ZXZ decomposition
% U = [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), A]*kron(H,eye(pow2(n-1)))*[B1,zeros(pow2(n-1));zeros(pow2(n-1)),B2]*kron(H, eye(pow2(n-1))) * [eye(pow2(n-1)), zeros(pow2(n-1));zeros(pow2(n-1)), C] ;
% Reference: Block-ZXZ
% http://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.22.034019
% We transform A = A2*A1', B1 = A1, B2 = A1*B, where A1, A2, and B are the
% original parameters in Block-ZXZ decomposition in above article.

    n = log(size(U,1))/log(2) ; 
    X = U(1:pow2(n-1),1:pow2(n-1)) ; 
    Y = U(1:pow2(n-1),1+pow2(n-1):end) ; 
    [UX,~,VX] = svd(X) ; 
    [UY,~,VY] = svd(Y) ; 
    WX = UX * VX'; WY = UY * VY'; 
    % transformed parameters A, B1, B2, C  
    C = -1j .* WX' * WY ; 
    B1 = X + Y*C' ; 
    B2 = X - Y*C' ;  
    A = ( U(1+pow2(n-1):end, 1:pow2(n-1)) + U(1+pow2(n-1):end, 1+pow2(n-1):end) * C' ) * B1' ; 
    
end
