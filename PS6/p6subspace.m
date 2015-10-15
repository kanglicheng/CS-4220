function [V, R] = p6subspace(A, m, maxiter, rtol)
% [V, R] = p6subspace(A, m, maxiter, rtol)
% Perform m-dimensional subspace iteration for matrix A,
% returning V (eigen matrix) and the upper triangular matrix R.
%
% The iterations continue till either the tolerance (rtol) has
% been satisfied or the maxiter iterations have been reached.
% 
% Convergence criterion is the Frobenius norm
% | AV_k - V_k R_k |_F  < rtol
%
n=length(A);
% pre-allocate V and X matrices (n by m):
V = zeros(n,m);
X = zeros(n,m);
%E = eyes(n,n);
% initialization: V_0:
[V, R]=qr(randn(n,m),0);
%V=eye(n,m);

k=1;
while (k < maxiter)
    X = A*V;
    % To get V at the next level k+1
    [V, R]=qr(X,0);
    % Schur's projection step
    Aproj=V'*(A*V);
    [U,T]=schur(Aproj);
    V=V*U;
    DIFF=A*V - X;
    residual=norm(DIFF,'fro');
    %residual=norm(DIFF);
    if rem(k,10)==0
        disp(sprintf('\nresidual=%g, k=%d', residual, k))
    end
    if ( residual < rtol )
        return
    end
    k=k+1;
end
end

