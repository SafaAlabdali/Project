%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using diagonal matrix DT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e,res]=Model_Jacobi_T(alpha , beta,N)
%create matrix 
Lx=1.0; % Physical size of the domain in X-direction
Ly=0.1; % Physical size of the domain in Y-direction
dx=Lx/N;
dy=Ly/N;
A=A2D(N,alpha,beta,dx,dy);
n=size(A,1);
% create random vector u (exact)
u=rand(n,1);
% find the vector b 
b=A*u;
%solve for Au=b
%u=jacobi(A,b);
% construct tridiagonal mat 
T = zeros(n,n);
T(1:1+n:n*n) = diag(A);
T(n+1:1+n:n*n) = diag(A,1);
T(2:1+n:n*n-n) = diag(A,-1);
u0=zeros(n,1);   % initial vector 

k=100; % number of iterations
res=zeros(k,1);  % residual 
e=zeros(k,1);    % error 
res_0=norm(A*u0-b);
res(1)=res_0;
e_0=norm(u0-u);
e(1)=e_0;

k=2;
while ( (res(k-1)/res_0) > 1E-4)
    r=b-A*u0;
    uk=u0+ T\r;        % jacobi(T,r);
    res(k)=norm(A*uk-b);
    e(k)=norm(uk-u); 
    u0=uk;
    k=k+1;
    if (k==100)
        break
    end
end

res=res/res_0;
e=e/e_0;
end 
