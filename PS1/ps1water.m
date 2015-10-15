function k=ps1water(omega, h)
%
% function k=ps1water(omega, h)
% finding k, the wave number, of the 
% shallow water waves dispersion equation
% using Newton's method.
% This function returns k (wave number, nondimensional) 
% by the given omega (frequency [1/s]) and h (depth [m]).
%
%
% TOL: tolerance for convergence: |F(k)| < TOL
% maxit: max number of Newton's iterations

TOL=1E-8;
maxit=100;
T_over_rho=7.2E-5;
g=9.8;

% anonymous functions defined  for f and f':
func = @(k) k*(g+T_over_rho*k^2)*tanh(k*h) - omega^2;

func_dot = @(k) (g+3*T_over_rho*k^2)*tanh(k*h) ...
  + k*(g + T_over_rho*k^2)*sech(k*h)^2* h;
%
  k = 0.1/h;  % initial guess based on kh << 1
  maxit = 100;  % maximum number of Newton's iterations

  for i=1:maxit
    if ( abs(func(k)) < TOL)
% we have got the root!
      return;
    else
      k = k - func(k)/func_dot(k);
      disp(sprintf('i=%d   x=%e    fx=%e', i, k, func(k)));
    end
  end
% something is wrong if we get here:
  disp('Something is wrong: quitting Newton"s method');
  k=NaN;
