%=====================================================
% This fuction uses SPAI algrithm  
%=====================================================

function M=SPAI(A,Apriori)
n=size(A,1);
I=speye(n);
M=sparse(n,n);

% loop that finds mk (columns of M)
for k=1:n 
  
    % find the non-zero rows of a priori pattern
    col_el=find(Apriori(k,:)); 
    A1=A(:,col_el);
    [rows_el,~]=find(A1);    % Find the position of non-zero element Per column wrt row
    rows_el=unique(rows_el); % Eliminate the repetitions
    n1=length(rows_el);
    n2=length(col_el);
    A1=A1(rows_el,:); %submatrix of size n1*n2
    %=====================================
    % compute the row indices I of the corresponding nonzero entries 
    %=====================================  
    e=I(:,n2);
    e=e(rows_el,:);
    %=====================================
    % compute th QR decomposition of A(I,J).
    %=====================================
    [Q R]=qr(A1);
    c=Q'*e;
    bj=R\c;
    %=====================================
    % construct the clumn of matrix M (mk)
    %=====================================
    mk=sparse(n,1);
    mk(col_el)=bj;
    M(:,k)=mk;
end

end
