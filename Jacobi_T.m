%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using diagonal matrix DT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,e,k]=Jacobi_T(A,b,N)

% construct tridiagonal mat 
T = zeros(N,N);
T(1:1+N:N*N) = diag(A);
T(N+1:1+N:N*N) = diag(A,1);
T(2:1+N:N*N-N) = diag(A,-1);

u0=zeros(N,1);   % initial vector 


res_0=norm(A*u0-b); % initial residual 
res = res_0;
k=1;
while ( (res(k)/res_0) > 1E-4)
    r=b-A*u0;
    uk=u0+ T\r;         % jacobi(T,r)
    res_k=norm(A*uk-b);
    u0=uk;
    res= [res res_k];   %vector that saves the residuals 
    k=k+1;
    if(k==2)
        break;
    end 
end
%%%%%%%%%%%%%%%%%%%%%%
u0=zeros(N,1);
% This part will be done to find the error 
u_ex=rand(N,1);
b_ex= A*u_ex;
e0= norm(u_ex -u0);

%%%%%%%%%%%%%%%%%

for i=1:k
    r =b_ex-A*u0;
    uk=u0+ T\r; %jacobi(T,r);
    e(i)=norm(uk-u_ex);
    u0=uk;
end
e =e/e0;
r=res/res_0;
end 