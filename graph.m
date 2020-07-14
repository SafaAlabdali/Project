%format shortEng
format compact

 alpha=1;
 beta=1E4;
 % residual 
 [r ,k] = Jacobi_D( alpha , beta);
 [r(1), r(10), r(50),r(100)]

 semilogy(1:k,r);       
  
 legend ('Residual');
 
 title(['alpha =',num2str(alpha),'Beta =',num2str(beta)])






