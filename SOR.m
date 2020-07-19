function [res , k ]= SOR(A,x0,b,w)


n=size(A,1);
res0=norm(A*x0-b);

res = res_0;
k=1;
while ( res(k)/res_0 > 1E-1)
    
 for i = 1:n
     sum= 0;
     for j=1:n
        if j~=i
        sum = sum + A(i,j)*x(j);
        end
     end
     x(i) = (1-w)*x(i) + w*(b(i)-sum)/A(i,i);
 end 
    resk=norm(A*x-b);
    res =[res resk];
    k=k+1;
    
end 

end 