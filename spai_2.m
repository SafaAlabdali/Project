function [M] = spai_2(A,S)
% Calculate the sparse approximate inverse M of A, for a given sparsity 
% pattern S. Here S is a matrix of the same size as A.

% Check that matrix sizes are valid and consistent
[m,n] = size(A);
assert(m==n,"A has to be a square matrix");
[m_s,n_s] = size(S);
assert((m_s==m)&&(n_s==n),"S has to have the same size as A");

% Check that A and S have no zero entries on diagonal
assert(all(logical(diag(A))),"A can not have zero diagonal entries");
assert(all(logical(diag(S))),"S can not have zero diagonal entries");

% Initialise M to zero matrix
M = sparse(m,n);
% Loop over columns of M
for j=1:n
    % find non-zero row-indices of column vector m_j with (m_j)_k = M_{kj}
    J_j = find(S(:,j));
    % find non-zero row-indices of all column vectors a_j with j in I_j,
    % where (a_j)_k = A_{kj}
    [rows,cols] = ind2sub(size(A(:,J_j)),find(A(:,J_j)));
    I_j = unique(rows);
    % Extract small submatrix A_j from A
    A_j = A(I_j,J_j);
    % Construct unit vector of the correct size
    e_j = zeros(size(I_j)); % initialise to zero
    e_j(find(J_j==j))=1;    % insert 1 at the correct position
    % Find min_{m_j} ||A_j.m_j - e_j||
    m_j = A_j\e_j;
    % Scatter elements back to matrix M, using the index set J_j
    M(J_j,j)=m_j;
end
end

