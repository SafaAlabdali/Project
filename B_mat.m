function B =B_mat(A)
n=size(A,1);
B=zeros(n,n);  % B=[B11 B12 ; B21 B22]

B(1:n/2,1:n/2)=A((n/2)+1:n,(n/2)+1:n); % B11
B(1:n/2,(n/2)+1:n)=A((n/2)+1:n,1:n/2); % B12
B((n/2)+1:n,1:n/2)= A(1:n/2,(n/2)+1:n);% B21
B((n/2)+1:n,(n/2)+1:n)= A(1:n/2,1:n/2);% B22

end 