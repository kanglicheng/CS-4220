function [y] = ps2mult(x)
%
% x is an (n by 1) vector
% A = reshape(1:n^2, n, n)
%
% This routine is to return y= A*x in O(n) flops
%
  n=length(x);
%
% the matrix  A = reshape(1:n^2, n, n);
% 
  a1=[1:n]'; 
%  w = ones(n,1);
% 
% Recognizing A (n by n) is a rank 2 matrix whose column 
% space is span{ a_1, w } (a1 and w are defined as above):
% 
% A = [a1  a1+n*w  a1+2n*w .... a1+(n-1)n*w ] 
%
% We exploit the structure of matrix A here so that  A*x 
% only uses  O(n) flops, instead of O(n^2) flops.
% 
%
%  sum1=0;
%  sum2=0;
%  for j=1:n
%   sum1=sum1+x(j);
%   sum2=sum2+(j-1)*x(j);
%  end
%
    sum1 = sum(x);    % the first term
    sum2 = (a1 -1)'*x; % the second term
    
% remember to multiply sum2 by n below (according to the simple formula
% we derived in the solution):
%
  y = sum1*a1 + n*sum2;