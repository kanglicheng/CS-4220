


function [ imf ]= p2lsqr( imblurd, H, niter)
%
%
% Deblur the input image using LSQR for n steps
% H: image transformation
% n: number of iterations
% 
% Output:
%  imf output image matrix
%
% This function uses the MATLAB built-in lsqr solver
% in solving the image deburring problem
% HV = V^{blur}.
%
%
TOL=1.0E-6;
%
% version 1: no preconditioning:
%

[m, n]=size(H);
%
% generate the column vectors (for three colors)
% from imblurd, and pre-allocate the solution vectors
%
for jj=1:3
    b{jj}=imblurd(:,jj);
    vv{jj}=zeros(n,1);
end
%
% call lsqr separately for each column vector.
%
for col=1:3
    vv{col}=lsqr(H, b{col}, TOL, niter);
end

imf=[vv{1} vv{2} vv{3}];
end