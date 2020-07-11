%format shortEng
format compact

% Create plots
t = tiledlayout(2,2); % Requires R2019b or later
title(t,'Beta =1')
for i=1:4
 k=10;
 alpha=1;
 beta=[0.01 ,0.1 ,1 ,10 ];
 [e,r] = Jacobi_schur( alpha , beta(i));
 [e(1), e(10), e(50),e(100)]
 [r(1), r(10), r(50),r(100)]
 
 
 nexttile
 semilogy(1:k,r);
 hold on 
 semilogy(1:k,e);
 legend ('Residual','Error');
 hold off
 title(['alpha =',num2str(beta(i))])
clc; clear all;
end 




