%read 
function A=read()
fileID = fopen('dumpMatrix000.ymat','r');
% Read record 1: 3 32bit integers
data = fread(fileID,3,'int32');
ns = data(1);        % system size 
na = data(2);         % number of diagonal entries 
nija = data(3);      %
nnz = nija - (na+1);  % number of ( non-zero off-diagonal )entries 
fprintf('na = %d nija = %d nnz = %d\n',na,nija,nnz);
% Read record 2: nija 32bit integers
rowPtr = fread(fileID,na+1,'int32');
colIdx = fread(fileID,nnz,'int32');
% Read record 3: nija-1 64bit doubles
values = fread(fileID,nija,'double');
% Read record 4: na 64bit doubles
rhs = fread(fileID,na,'double');
fclose(fileID);
A = diag(values(1:na));
for i=1:na
for k=rowPtr(i):rowPtr(i+1)-1
A(i,colIdx(k-na)+1) = values(k+1);
end
end
%plot(1:size(rhs),rhs)
spy (A)

end 