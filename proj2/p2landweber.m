function [ imresult ] = p2landweber(imblurd, H, niter, alpha)
%
% imblurd: blurred image (N by 3)
% H: image transformation matrix (N by N)
% niter: number of fixed point iterations
% alpha: the Landweber iteration parameter (a small parameter)
%
% Using Landweber interations for solving the image
% deburring problem
[m,n]=size(H);
AHT_Vb=alpha*H'*imblurd;
MM=speye(n)-alpha*H'*H;
% 
VV_Landweber=zeros(n,3);
%
% VV_Landweber stores the matrix current (k)
% and updated  level (k+1):
%
for jj=1:niter
   VV_Landweber=MM*VV_Landweber+AHT_Vb;
end
imresult=VV_Landweber;
end