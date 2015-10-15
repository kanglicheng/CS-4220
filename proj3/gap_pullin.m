function [u,V] = gap_pullin(pgap)
% [u,V] = gap_pullin(pgap)
%
% Compute the pull-in displacement and voltage for the given
% structure.
%
% Inputs:
%    pgap: Structure returned by gap_setup
%
% Outputs:
%    u: Displacement vector at pull-in state
%    V: Pull-in voltage
%
% === Load gap finite element code

%pgap = gap_setup('c', 80);

% nelt = pgap.nelt;
% Ce = pgap.Ce;
% M = pgap.M;
% K = pgap.K;
% N = pgap.N;
% wg = pgap.wg;
% Itip = pgap.Itip;

ndof = pgap.ndof;

% the following interval [disp1, disp2] must
% bracket the transition point of the curve
%
disp1=-0.8*pgap.g;
disp2=0;

%
%
V1=10;
V2=0.01;

% use Bi-Section method to find the pull-in voltage/displacement
% by using the test that the Jacobian [K -Ju] 
% changes from positive definite (p==0)
% to indefinite (p>0). We change p==0 to p<0 for finding this root:
% we start from a=disp1 (lower bound) and b=disp2 (upper bound).
% | V_pullin - p | < atol
% Bisection code below:
%
a=disp1;
b=disp2;
u=zeros(ndof,1);
u2=zeros(ndof,1);
fa=0; fb=0;
[fa, V1, u]=gap_matrix_PDtest(pgap, u, V1, a);
[fb, V2, u2]=gap_matrix_PDtest(pgap, u2, V2, b);
if (a >= b) | (fa*fb > 0)
    fprintf('problem with the initial input:quitting\n');
    V=NaN;
    %u=zeros(ndof,1);
    return;
end
atol=1e-5;
n = ceil( log2(disp2-disp1) - log2(2*atol));
V=V1;
for k=1:n+1
   p=(a+b)/2; 
   % p is the next displacement value for bisection method:
   [fp, V, u]=gap_matrix_PDtest(pgap, u, V, p);
   if (fa*fp < 0)
       b=p;
       fb=fp;
   else
       a=p;
       fa=fp;
   end
end
%
% end of bisection
%
%fprintf('\n***  The pull-in voltage is V=%g', V);
%fprintf('\n***  using atol=%g and bisection method\n', atol);

end
