% This fuction solves linear system Av=r 
% by Schur complement (basic)
function v=solve (A,r)
n=size(A);
n=n(1);
r1=r(1:n/2);
r2=r((n/2)+1:n);

%================ case of matrix B =================
%A=B_mat(A);


%=================extracting A to four blocks =================
A11=A(1:n/2,1:n/2);
A12=A(1:n/2,(n/2)+1:n);
A21=A((n/2)+1:n,1:n/2);
A22=A((n/2)+1:n,(n/2)+1:n);

%==================choice of A22_tild========================
m=size(A22,1);
A22_tild= zeros(m,m);
%= 1) A22 itself
%A22_tild=A22;

%= 2) diagonal of A22
%A22_tild=diag(diag(A22)); 

%= 3) identity
%A22_tild=eye(size(m));   

%= 4) tridiagonal
%A22_tild(1:1+m:m*m) = diag(A22);   
%A22_tild(m+1:1+m:m*m) = diag(A22,1);
%A22_tild(2:1+m:m*m-m) = diag(A22,-1);
%===========================================
% inv(A22_tilde)
M=sappinv(A22,Apriori(A22));

%=================compute r1_tild==========
r0 = M* r2; %jacobi(A22_tild,r2);
% r0 = A22_tild\r2; (use this for testing only)
r1_tild=r1- A12 * r0;

%================solve sv=r1_tild=============
S = A11 - A12 * M *A21; %inv(A22_tild) *A21 ;
% jacobi method 
% v1 = S\r1_tild; (use this for testing only)
v1=jacobi(S,r1_tild);
%================compute v2_tild=======
v2=jacobi(A22_tild,(r2-A21*v1));
%v2 = A22_tild\(r2-A21*v1); (use this for testing only)

v=zeros(n,1);
v(1:n/2) = v1;
v(n/2+1:n) = v2 ;


end 