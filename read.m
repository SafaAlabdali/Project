%read 
function [A,rhs,na]=read()
matfile = fopen('dumpMatrix000.ymat','r');
% Read record 1: 3 32bit integers
data = fread(matfile,3,'int32');
ns = data(1);         % system size 
na = data(2);         % number of diagonal entries 
nija = data(3);       % size of matrix values val and row pointers/column indices ija
nnz = nija - (na+1);  % number of ( non-zero off-diagonal )entries 

% Read record 2: nija 32bit integers
rowPtr = fread(matfile,na+1,'int32');%Row pointers, length=number of rows (ie. number of diagonal entries) + 1
colIdx = fread(matfile,nnz,'int32'); %Column indices, length=total number of non-zero entries=na+nnz
% Read record 3: nija-1 64bit doubles
values = fread(matfile,nija,'double');
% Read record 4: na 64bit doubles
rhs = fread(matfile,na,'double');
fclose(matfile);
A = sparse(na,na);
for k=1:na
    A(k,k) = values(k);
end
for i=1:na
for k=rowPtr(i):rowPtr(i+1)-1
A(i,colIdx(k-na)+1) = values(k+1);
end
end
%plot(1:size(rhs),rhs)
%spy (A)

end 