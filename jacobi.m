%This function solves linear system Av=r by 
% jacobi method 

function v= jacobi(A,r)
n=size(A);
n=n(1);
D=diag(diag(A)); % diagonal matrix 
v=zeros(n,1); % initial vector 
k=2; % number of iteration  
for i=1:k
    vk = v +( D \(r-A*v) );
    v=vk;
end 


end 