%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using diagonal matrix DT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,k]=Jacobi_D(alpha,beta)

[A,b,N]=read(); % read matrix and RHS 
Lx=1.0;         % Physical size of the domain in X-direction
Ly=0.1;         % Physical size of the domain in Y-direction
dx=Lx/N;
dy=Ly/N;

% construct tridiagonal mat 
T = zeros(N,N);
T(1:1+N:N*N) = diag(A);
T(N+1:1+N:N*N) = diag(A,1);
T(2:1+N:N*N-N) = diag(A,-1);

u0=zeros(N,1);   % initial vector 
  

res_0=norm(A*u0-b); % initial residual 
res = res_0;
k=1;
while ( res(k) > 1E-4)
    r=b-A*u0;
    uk=u0+ T\r;         % jacobi(T,r)
    res_k=norm(A*uk-b);
    u0=uk;
    res= [res res_k];   %vector that saves the residuals 
    k=k+1;
end
res=res/res_0;
end 