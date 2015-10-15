function v = inverseIteration(A, lambda, tolerance)
%v = inverseIteration(A, lambda, tolerance)
sigma=lambda*(1. - 1E-5);
maxiter=20;
k=1;
n=length(A);
v=rand(n,1); % v_0 has to be non-zero
v_new = zeros(n,1); % pre-allocate v_new

AA = A - sigma*eye(n,n);
[L, U]=lu(AA);

while ( k < maxiter )
  % we need AA^{-1} * v:
  v_new = U\(L\v);
  v_new = v_new/norm(v_new);
  resid=norm(v_new - v);
  if ( resid < tolerance )
    v=v_new;
    return
  end
  disp(sprintf('k=%d, residual=%e', k, resid))
  k=k+1;
  v=v_new;
end
end

