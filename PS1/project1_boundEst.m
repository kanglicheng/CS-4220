%Generating a test sparse matrix for CS4220 Project1
%
%Prob=UFget(2317);
%A=Prob.A;
%
%Load the matrix from roadNet-CA.mat file
%
%
Z=load('roadNet-CA.mat');
A = Z.Problem.A;
alpha = 0.9;

%testing one (s,t) pair
s=180; t=1;

n = length(A);
d = full(sum(A))';
INZ = find(d);
n = length(INZ);
d = d(INZ);  % get rid of the zero-degree nodes
A = A(INZ, INZ);
D = spdiags( d, 0, n, n );
Dinv = spdiags( 1./d, 0, n, n);
N = D -alpha*A;

T = A*Dinv;
M=speye(n) - alpha*T;
p=2;  %This problem is a rank-2 update to matrix M

%%%%%%%%%%%%%%%%%%%%%%%

global yy I_p;  % yy is used repeatedly in Sherman_Morrison_Woodbury()
I = speye(n);
I_p=speye(p);
e{s}=I(:,s); 
e{t}=I(:,t);
%
% R'*R = S'*N*S, so Nw=b can be solved as
%        w = S*(R\(R'\(S'*b)));
%
[R, pp, S]=chol(N);
%ww=S*(R\(R'\(S'*e{t}))); % this is an error! should be e{s}
ww=S*(R\(R'\(S'*e{t})));
% PMQ=tL*tU
%[tL, tU, P, Q]=lu(M);
%yy=Q*(tU\(tL\(P*e{s})));  

%A is symmetric, take its upper triangular part for (a,b)
[r,c]=find(triu(A));

Nedges=length(r);
fab_bound=zeros(Nedges,1); 
temp1=1/d(t)*(1+alpha)/(1-alpha);
%deg_two = find(d>1);
%Nedges2=length(deg_two);
%fab_bound(deg_two)=temp1*( ww(r(deg_two))./(d(r(deg_two))-1)+...
%                   ww(c(deg_two))./(d(c(deg_two))-1) );
ind=(1:Nedges)';
fab_bound(ind)=temp1*( ww(r(ind))./(d(r(ind))-1)+...
                  ww(c(ind))./(d(c(ind))-1) );
% It is critically important to get rid of Inf and NaN from
% the fab_bound (since we have not filtered out single degree nodes:

fab_bound( isinf(fab_bound) | isnan(fab_bound) )=0;

[opt, indsort] = sort( fab_bound, 'descend' );
[optedge, indexx] = max( fab_bound );
a=r(indsort(1)); % a, b are the optimal edge pair
b=c(indsort(1)); %
disp(sprintf('max. bound occurs at %d, a=%d, b=%d', indsort(1), a, b))
%disp(sprintf('\n (a=%d,  b=%d)',  a, b))

% e{a}, e{b}, u{1} etc. are cell arrays, which
% can be used for arrays of vectors or matrices.
% This is equivalent to array of any data structure
% in Java, C++ and C.

e{a}=I(:,a);
e{b}=I(:,b);

if ( d(a) > 1 )
    u{1}=alpha/d(a)*A(:,a) - (A(:,a) - e{b})*alpha/(d(a)-1);
else
    u{1}=alpha/d(a)*A(:,a);
end

if ( d(b) > 1 )
    u{2}=alpha/d(b)*A(:,b) - (A(:,b) - e{a})*alpha/(d(b)-1);
else
    u{2}=alpha/d(b)*A(:,b);
end

U = [u{1} u{2}];
V = [e{a} e{b}];  % V' is the 2xn matrix
%
% We need not generate M_hat in this code. Sherman-Morrison-
% Woodbury produces M_hat^{-1} for us.
%
%  M_hat = M+U*V';
% (M + U*V')^{-1} is obtained by the following:
%
% PMQ=tL*tU
[tL, tU, P, Q]=lu(M);
yy=Q*(tU\(tL\(P*e{s})));

x=Sherman_Morrison_Woodbury(tL, tU, P, Q, U, V');
f=e{t}'*x;

fval=full(f);
disp(sprintf('\n for t=%d,  s=%d ', t, s));
disp(sprintf('the optimal edge is (%d  %d),   f=%g', a, b, fval));
     
 
    