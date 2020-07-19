function C = Apriori (A)

n=size(A,1);
C=zeros(n,n);

for i=1:n     % for each raw
    
    for j=1:n % at row i and jth entry
       if A(i,j)~= 0
           C(i,j)= 1;
       else 
           C(i,j)= 0;
       end
    end   
end 
end 