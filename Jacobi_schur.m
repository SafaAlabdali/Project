%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program test Jacobi Method using tridiagonal matrix T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res ,e ,k]=Jacobi_schur(A,b)
[N,N_]= size(A);
%=================extracting A to four blocks =================
n=size(A,1);
A11=A(1:n/2,1:n/2);
A12=A(1:n/2,(n/2)+1:n);
A21=A((n/2)+1:n,1:n/2);
A22=A((n/2)+1:n,(n/2)+1:n);
%==============================================================
u0=zeros(N,1);      % initial vector 
res_0=norm(A*u0-b); % initial residual 
res = res_0;
k=1;
while ( (res(k)/res_0) > 1E-6)
    r=b-A*u0;
   
    uk=u0+ solve (A11,A12,A21,A22,r,n);%solve(A,r);
    
    res_k=norm(A*uk-b);
    u0=uk;
    res = [res res_k];   %vector that saves the residuals 
    k=k+1;
    if(k==5)
       break;
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%
u_ex=rand(N,1);
b_ex=A*u_ex;
u0=zeros(N,1);
e0=norm(b_ex-A*u0);

for i=1:k
    r=b_ex-A*u0;
    uk=u0+ solve (A11,A12,A21,A22,r,n);
    e(i)=norm(uk-u_ex);
    u0=uk;
end 

e=e/e0;
res=res/res_0;

end 
