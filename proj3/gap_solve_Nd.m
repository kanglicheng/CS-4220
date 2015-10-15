function [u, V, rnorms] = gap_solve_Nd(pgap, u, V, d, rtol, maxiter)
% [u, V, rnorms] = gap_solve_Nd(pgap, u, V, d, rtol, maxiter)
%
% Attempt to solve the gap equation by Newton's method.  Applies
% displacement control to the tip.
%
% Inputs:
%   pgap: Structure returned by gap_setup
%   u: Initial guess of displacement vector
%   V: Initial guess of voltage
%   d: Tip displacement
%   rtol: Residual tolerance for convergence (default: 1e-8)
%   maxiter: Maximum iterations allowed (default: 10)
%
% Outputs:
%   u: Computed displacement vector
%   V: Computed voltage level
%   rnorms: Residual norms at each iteration
%

  % Assign default tolerances if not specified
  if nargin < 5, rtol = 5e-7;  end
  if nargin < 6, maxiter = 15; end
  
  %pgap=gap_setup('c',80);
  
  %maxiter=10;
  % initial V
  %V=1;
  %d=-4.0464e-3;
  %d=-0.33022
  
  K=pgap.K;
  Itip=pgap.Itip;
  ndof=pgap.ndof;
  u=zeros(ndof,1);
  
  %[u, rnorms1] = gap_solve_NV(pgap, u, V, rtol, 15);

  %set etip, the unit base vector for the tip displacement
  etip=zeros(ndof,1);
  etip(Itip)=1;  
  
  nzero=zeros(ndof-1,1);
  
  % Newton's method is applied to the 
  % vector equations as follows:
  %
  % Fr=[ Kn*z - f_mod ] =0,
  % where z = [u_1 u_2 . . .   u_(Itip-1) u_n V^2]' ,  
  % with u(Itip) being removed from u, and V^2 being
  % added to the tail.
  % The critical step is 
  % to modify the RHS vector as follows:
  % Take out the entire column of 
  % K(:, Itip) and move d*K(:,Itip) to the right:
  %  so f_mod = f(u, V) - d*K(:,Itip) 
  %
  F=zeros(ndof,1);
  % Kn is obtained by taking out the Itip column K(:,Itip),
  % and fill in the last column by zeros (for V^2 is the
  % new variable added to the vector, but V^2 has no
  % explicit contribution to the equations:
  Kn=[K(:,1:Itip-1) K(:,ndof) zeros(ndof,1)];
  % Kr is by taking K(:,Itip) out from K 
  Kr=[K(:,1:Itip-1) K(:,ndof)];
  rnorms=zeros(maxiter,1);
  % Note that u is still ndof by 1,
  % u(Itip) is out, V^2 is added to the tail.
  % Note that we do need the original u whenever
  % we evaluate fe(u, V).
  % z is the new solution vector as follows:
  % z = [u_1 u_2 . . .   u_(Itip-1) u_n V^2]' ,
  %
  z=[u(1:Itip-1);u(ndof);V^2];
  %fprintf('\n*** initial V=%e, d=%e\n', V, d);
  % Main iteration ***
  for iter=1:maxiter
    % call gap_force() for f, and Jacobian matrices Ju and JV.
    % obtain V from the solution vector z:
    V=sqrt(z(end));
    % Be careful to reconstruct u correctly:
    u=[z(1:Itip-1);d;z(ndof-1)];
    [f, Ju, JV]=gap_force(pgap, u, V);
    % f_mod is the right-hand side vector after
    % we move u(Itip)=d over
    f_mod=f-d*K(:,Itip);
    % ur=[u(1:Itip-1); u(end)];
    % form the F function for Newton:
    F=Kn*z-f_mod; 
    % calculate the residual norm of F:
    rnorms(iter)=norm(F);
    % is it converged yet?
    if rnorms(iter) < rtol
        rnorms=rnorms(1:iter);  %resize rnorms for output
        return;
    end
    % solve J_F*p_k = -F
    % Construct the  Jacobian matrix:
    Jur=[Ju(:,1:Itip-1) Ju(:,Itip+1)];
    %JF=[Kn(:,1:ndof-1)-Jur -JV];
    JF=[Kr-Jur -1*JV];
    %pvec=-pinv(JF)*F;
    pvec=-JF\F;
    z=z+pvec;
    %fprintf('\niter=%d, rnorms=%e', iter, rnorms(iter));
    %fprintf('\ntip Disp= %e, V=%e\n', u(Itip), sqrt(z(end)));
  end
  %
  % If we got here, we didn't converge in time
  warning('Newton iteration did not converge')
  warning('iter=%d, residual norm =%e, V=%e', iter, rnorms(iter), V)
end