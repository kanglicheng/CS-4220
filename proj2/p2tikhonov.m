function imresult = p2tikhonov(imblurd, H, beta)
% Solving the Tikhonov Reguarlization problem
% input: 
%  imblurd: matrix (image) (n by 3)
%  H: the image transform matrix (n by n) 
%  beta: the regularization parameter
% output:
%  imresult: output image (n by 3)
%  
%  Solving an ill-conditioned least-squares problem
%  by using Tikhonov regularization
%    V^{tik} = argmin_V |HV -Vb|^2 + beta^2*|V|^2
%  where Vb is the blurred image, H is the image
%  transformation matrix, V is the solution, and beta
%  is the regularization parameter for solving this 
%  kind of ill-conditioned problem.

[m, n]=size(H);
N=H'*H+beta^2*speye(n);
%[L,U,P,Q]=lu(N);
%imresult=Q*(U\(L\(P*H'*imblurd)));
%
% sparse Cholesky decomposition
% S is the permutation matrix
%
[R, p, S]=chol(N);
imresult=S*(R\(R'\(S'*H'*imblurd)));
end

