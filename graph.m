%format shortEng
format compact

 tic
 [A,b,N]=read();
 [r ,k] = Jacobi_schur(A,b,N)  % it returns residual and number of iterations to reach specific tolerance 
 toc
 [r(1), r(10), r(50),r(100)]

 semilogy(1:k,r);       
  
 legend ('Residual');
 
 title(['The residual up tp ',num2str(k),'iterations'])






