% Source:
%   http://www.itl.nist.gov/div898/strd/nls/data/thurber.shtml
%
function ps8thuber
xy =[80.574E0      -3.067E0
    84.248E0      -2.981E0
    87.264E0      -2.921E0
    87.195E0      -2.912E0
    89.076E0      -2.840E0
    89.608E0      -2.797E0
    89.868E0      -2.702E0
    90.101E0      -2.699E0
    92.405E0      -2.633E0
    95.854E0      -2.481E0
    100.696E0      -2.363E0
    101.060E0      -2.322E0
    401.672E0      -1.501E0
    390.724E0      -1.460E0
    567.534E0      -1.274E0
    635.316E0      -1.212E0
    733.054E0      -1.100E0
    759.087E0      -1.046E0
    894.206E0      -0.915E0
    990.785E0      -0.714E0
    1090.109E0      -0.566E0
    1080.914E0      -0.545E0
    1122.643E0      -0.400E0
    1178.351E0      -0.309E0
    1260.531E0      -0.109E0
    1273.514E0      -0.103E0
    1288.339E0       0.010E0
    1327.543E0       0.119E0
    1353.863E0       0.377E0
    1414.509E0       0.790E0
    1425.208E0       0.963E0
    1421.384E0       1.006E0
    1442.962E0       1.115E0
    1464.350E0       1.572E0
    1468.705E0       1.841E0
    1447.894E0       2.047E0
    1457.628E0       2.200E0];

% Data points
y = xy(:,1);
x = xy(:,2);

% Mesh for plotting f
xx =linspace(x(1), x(end));

% Initial values
b = [1000;
    1000;
    400;
    40;
    0.7;
    0.3;
    0.03];

% Certified values
bref = [1.2881396800E+03;
    1.4910792535E+03;
    5.8323836877E+02;
    7.5416644291E+01;
    9.6629502864E-01;
    3.9797285797E-01;
    4.9727297349E-02];

% Tolerance on norm(J'*r)
rtol = 1e-8;

% Record residual norms
resids = [];

% Evaluate initial residual
n = b(1) + (b(2) + (b(3) + b(4)*x).*x).*x;
d =    1 + (b(5) + (b(6) + b(7)*x).*x).*x;
f = n./d;
r = f-y;

% TODO: Iterate to find b s.t. norm(r)^2 is minimal
%       I recommend Gauss-Newton with line search;
%       you may prefer Levenberg-Marquardt or something fancier.
maxiter=80;
ncol=length(b); % b is the unknown vector of beta coefficient
fp={};  % fp is a cell array: array of vectors in the Jacobian matrix

%evaluate the initial Jacobian matrix
fp{1}= 1 ./ d;
fp{2}= fp{1}.*x;   %fp{2}=x ./ d;
fp{3}= fp{2}.*x;   %fp{3}= x.*x ./ d;
fp{4}= fp{3}.*x;   %x.^3 ./ d;

fp{5}=-n ./ d.^2 .*x;
fp{6}=fp{5}.*x;
fp{7}=fp{6}.*x;

J = [ fp{1} fp{2} fp{3} fp{4} fp{5} fp{6} fp{7} ];
alpha=0.5; % line search will provide the value of  alpha
% main Gauss-Newton iteration
for iter=1:maxiter
    % J and r are either from the intialization or updated after
    % solution (b) is updated at the (k+1) level later:
    %
    % checking for convergence here:
    resids(iter)=norm( J'*r );
    if (resids(iter) < rtol )
        break;
    end
    
    % Solving Gauss-Newton equation for search direction p:
    p = (J'*J)\(-J'*r);
    
    [phi]=phi_beta(b, x, y);
    
    grphi=J'*r;
    
    % solution b is updated via a line search:
    [bn,alpha]=lineSearch(b, p, x, y, phi, grphi);
    b=bn;
    fprintf('\n iter=%d, alpha=%g, res=%g\n', iter, alpha, resids(iter));
    %b = b + alpha*p;
    % update n, d and r with the new b (beta):
    n = b(1) + (b(2) + (b(3) + b(4)*x).*x).*x;
    d =    1 + (b(5) + (b(6) + b(7)*x).*x).*x;
    r = n./d - y;
    %
    % calculate the updated columns of the Jacobian matrix:
    % fp{i}=\partial f_i / \partial\beta:
    % fp{} is a cell array data:
    %
    % vectorization by the term-by-term operators
    fp{1}= 1 ./ d;
    fp{2}= fp{1}.*x;   %fp{2}=x ./ d;
    fp{3}= fp{2}.*x;   %fp{3}= x.*x ./ d;
    fp{4}= fp{3}.*x;   %x.^3 ./ d;
    
    fp{5}=-n ./ d.^2 .*x;
    fp{6}=fp{5}.*x;
    fp{7}=fp{6}.*x;
    
    %J = [ fp{1} fp{2} fp{3} fp{4} fp{5} fp{6} fp{7} ];
    for col=1:ncol
        J(:,col)=fp{col};
    end
    
end

% Check consistency with reference values
fprintf('Relative error vs ref values\n');
fprintf('%d: %e\n', [1:7; abs((bref-b)./bref)']);

% Plot results
figure(1);
nxx = b(1) + (b(2) + (b(3) + b(4)*xx).*xx).*xx;
dxx =    1 + (b(5) + (b(6) + b(7)*xx).*xx).*xx;
fxx = nxx./dxx;
plot(x, y, '.', xx, fxx);

% Plot convergence
figure(2);
semilogy(resids);
%
%
end
% **** Additional functions for this code ****
function [bn,alpha]=lineSearch(b, p, x, y, phi, grphi)
% [bn]=lineSearch(b, bnew, p, J, f, r)
%   A "weak line search" method for Gauss-Newton in nonlinear LSQ.
%
%   inputs:
%           b: beta at iteration (k)
%           bnew : beta at iteration (k+1)
%           p: Gauss-Newton search direction (descent)
%           phi: phi at k-level
%           grphi:  grad(phi) at k-level
%   output:
%           bn: beta at (k+1) level
%   (Taken from Ascher and Greif: page 264)
%
alphamax=0.8; 
alphamin=0.01; 
sigma=1E-3;

grpsi=p'*grphi;
alpha=alphamax;
%initially try max alpha
bn=b+alpha*p; 
[phi_n]=phi_beta(bn,x,y);

while (phi_n > phi+sigma*alpha*grpsi)*(alpha > alphamin)
    mu=-0.5*grpsi*alpha/(phi_n-phi-alpha*grpsi);
    if mu < 0.1
        mu=0.5;
    end
    %try this new alpha
    alpha=mu*alpha;
    bn=b+alpha*p;
    [phi_n]=phi_beta(bn,x,y);
end
end
