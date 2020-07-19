%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using diagonal matrix DT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,k]=Jacobi_D(A,b,N)

D=diag(diag(A)); % diagonal matrix 
u0=zeros(N,1);   % initial vector 
  
res_0=norm(A*u0-b); % initial residual 
res = res_0;
k=1;
while ( (res(k)/res_0) > 1E-1)
    r=b-A*u0;
    uk=u0+ D\r;         % jacobi(D,r)
    res_k=norm(A*uk-b);
    u0=uk;
    res= [res res_k];   %vector that saves the residuals 
    k=k+1;
    if(k==100)
        break;
    end 
end
res=res/res_0;
end 