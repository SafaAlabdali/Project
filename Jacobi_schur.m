%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using tridiagonal matrix T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res ,k]=Jacobi_schur(alpha,beta)

[A,b,N]=read(); % read matrix and RHS 
Lx=1.0;         % Physical size of the domain in X-direction
Ly=0.1;         % Physical size of the domain in Y-direction
dx=Lx/N;
dy=Ly/N;

u0=zeros(N,1);   % initial vector 
res_0=norm(A*u0-b); % initial residual 
res = res_0;
k=1;
while ( res(k) > 1E-4)
    r=b-A*u0;
    uk=u0+ solve(A,r);       
    res_k=norm(A*uk-b);
    u0=uk;
    res= [res res_k];   %vector that saves the residuals 
    k=k+1;
end

res=res/res_0;

end 
