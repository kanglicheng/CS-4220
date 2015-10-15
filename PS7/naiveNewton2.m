% PS7 Problem 2
% rtol=1E-12;
  n=2;
  rtol=1E-12;
  kmax=20;
  Fnorm=zeros(kmax,1);
  x=ones(n,1);  % initialization x0=(1,1)
  for k = 1:kmax
    % x is x(1), y is x(2) in the following expressions:
    f1=x(1)^2+x(1)*x(2)^2-9;
    f2=3*x(1)^2*x(2)-x(2)^3-4;
    F = [f1 f2]';
    Fnorm(k)=norm(F);
    if norm(F) < rtol
      break;
    end
    J = [2*x(1)+x(2)^2  2*x(1)*x(2);  6*x(1)*x(2)  3*x(1)^2-3*x(2)^2];
    dx = J\F;
    x = x-dx;
    disp(sprintf('\n k=%d,  x=[%e  %e]',k,x(1),x(2)))
    disp(sprintf('\n normF = %e', norm(F)));
  end
  
  hx=[1:k];
  semilogy(hx, Fnorm(1:k));
  
  