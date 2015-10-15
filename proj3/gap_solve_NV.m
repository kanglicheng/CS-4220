function [u, rnorms] = gap_solve_NV(pgap, u, V, rtol, maxiter)
% [u, rnorms] = gap_solve_NV(pgap, u, V, rtol, maxiter)
%
% Attempt to solve the gap equation by Newton's method.
%
% Inputs:
%   pgap: Structure returned by gap_setup
%   u: Initial guess of displacement vector
%   V: Voltage
%   rtol: Residual tolerance for convergence (default: 1e-8)
%   maxiter: Maximum iterations allowed (default: 10)
%
% Outputs:
%   u: Computed displacement vector
%   rnorms: Residual norms at each iteration
%
%
  % Assign default tolerances if not specified
  if nargin < 4, rtol = 1e-8;  end
  if nargin < 5, maxiter = 10; end
  %
  %
  K=pgap.K;
  %
  % Main iteration of Newton's method for F(u,V)=0:
  %   F(u,V)=Ku -f_e =0
  %
  F=zeros(pgap.ndof,1);
  rnorms=zeros(maxiter,1);
  for iter=1:maxiter
    % call gap_force() for f, and Jacobian Ju:
    [f, Ju]=gap_force(pgap, u, V);
    % calculate the residual norm of F=K*u-f
    F=K*u-f;
    rnorms(iter)=norm(F);
    % is it converged yet?
    if rnorms(iter) < rtol
        rnorms=rnorms(1:iter);  %resize rnorms for output
        return;
    end
    % solve J*p_k = -F = -(K*u_k - f)
    JF=K-Ju;
    pvec=-JF\F;
    u=u+pvec;
  end
  
  % If we got here, we didn't converge in time
  warning('Newton iteration did not converge')
  warning('iter=%d, residual norm =%e, V=%e', iter, rnorms(iter), V)
  
end