%format shortEng
format compact

 tic
 [A,b,N]=read();
 %============= Schur complement ===========
 %A=B_mat(A);
 %A = sparse(A);
 %====================
 [r,e,k] =Jacobi_schur(A,b);  % it returns residual and number of iterations to reach specific tolerance 
 toc
% [r(1), r(10), r(50),r(100)] 
tiledlayout(2,1)
nexttile
semilogy(2:k,r(2:k));
title('Residual')
nexttile
semilogy(2:k,e(2:k));
title('Error')
title(['The residual and Error up tp ',num2str(k),' iterations']);