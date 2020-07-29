%format shortEng
format compact

% Create plots
t = tiledlayout(2,2); % Requires R2019b or later
title(t,'Depened on experiment')

for i=1:4
 N=32;
alpha=1;
 beta=[0.0001 ,0.001,10 ,1000 ];
 
 [e,r] = Model_Jacobi_T(alpha , beta(i),N);
 
 % extracting the number of iteration 
 k=size(r); 
 k=k(1);
 %[e(1), e(10), e(50),e(100)]
 %[r(1), r(10), r(50),r(100)]
 nexttile
 semilogy(1:k,r);
 hold on 
 semilogy(1:k,e);
 legend ('Residual','Error');
 hold off
title(['beta =',num2str(beta(i))]);
clc; clear all;
end 




