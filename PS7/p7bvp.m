% for Problem 3 of pset7:
% Solve:
%   v" + \gamma exp( v ) =0, 0< x <1,
%    v(0)=v(1)=0
% (by the finite-difference method)
% n equally spaced intervals of length h
n=100;
gmax=100;
h=1/(n+1);
tol=1e-10; % tolerance for controling the Newton's iterations
% don't get confused: x is not the solution
% vector, x is the node vector in (0,1), 
% v(x) is the solution vector
x=0:h:1; 
gammas=linspace(1,3.5,gmax);
e=ones(n-1,1);
% T is the Tridiagonal matrix of n by n from
% which we can build J(v) and F(v):
%
T=diag(e,-1)+diag(e,1)-2*eye(n); 
v=zeros(n,1);  %v is the solution vector!
% the following vv matrix is for plotting the solution v 
% for different gammas:
% n+2 because solution vector plus two nodes: v(0) and v(1)
vv=zeros(n+2, gmax);  % vv is a matrix for plotting all solutions
lmax=zeros(gmax,1);
maxiter=30;
k=1;

for k=1:gmax
    for j=1:maxiter
        bb=gammas(k)*exp(v);
        % The FDM approximation gives the following 
        % equation F=0. (F is an n by 1 vector, note that 
        %                the mesh has 102 nodes.)
        F=(1/h^2)*T*v + bb ;
        J=(1/h^2)*T + diag(bb);
        %
        % Newton's method is to solve
        % F + J*p = 0, so 
        % p= -J\F
        %
        p = -J\F;
        v = v + p;
        norm_p = norm(p);
        norm_v = norm(v);
        if ( norm_p < tol*(1+norm_v) )
            fprintf('\n***j=%d,   norm_v=%e\n', j, norm_v)
            %largest eigenvalue of J(x*): eig which is closest to zero
            %because J is negatively definite, so use 'sm' option
            %in eigs function. J is approaching singular for eig is 
            %rapidly approaching zero.
            lmax(k)=eigs(J,1,'sm'); 
            vv(:,k)=[0 v' 0]';
            break;
        end
    end
end
% Plot the solution from the vv matrix versus x:
figure(1);
plot(x, vv);
%
% Plot the lambda_max of J(v*)
%
%figure(2);
%plot(gammas,lmax);
