function [x] = Sherman_Morrison_Woodbury(tL, tU, P, Q, U, V)
% Given a sparse LU decomposition of matrix M
% [tL, tU, P, Q]=lu(M);
% U (n by p) and V (p by n) are the rank-p update matrices
% such that M_hat = (M + UV), a rank-p perturbation of M 
% Use Sherman-Morrison-Woodbury to find 
% x=(M + PV)^{-1} b
% For efficiency purposes,
% yy, a global vector, 
%  M^{-1} yy = b,  yy = Q*(tU\(tL\(P*b)))
% has been precomputed once and used many times.
% 
% This function returns x=(M+UV)^{-1}*b
%
global yy I_p;
% To improve efficiency, p no longer
% generated within the function. I_(pxp) is
% a global variable.
% [n,p]=size(U);
%  yy is declared as global and precomputed before calling
%     this function
%%% yy=Q*(tU\(tL\(P*b)));
%
W=Q*(tU\(tL\(P*U)));
%
%C=speye(p)+V*W;
%
C=I_p+V*W;
z=C\(V*yy);
x=yy-W*z;
end

