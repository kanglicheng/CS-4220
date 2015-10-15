%
% Testing script for the column-stochastic matrix A:
%
N=50000;
P=[0 1/2 1/2 1/2 0 0;
   1/2 0 1/2 0   0 0;
   1/2 1/2 0  0  0 0;
   0  0   0   0  0 0;
   0  0   0   1/2  0 1;
   0  0   0   0  1 0  ]
[n,n]=size(P);
e=ones(n,1);
alpha=0.85;
A=alpha*P +(1-alpha)/n*e*e';

%v=e;
%for i=1:N
%   v = A*v; 
%end