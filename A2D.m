function M= A2D(N,alpha,beta,dx,dy)

% -\Laplace u + \alpha u + \beta (u-v)=r_u
% -\Laplace v + \alpha v + \beta (v-u)=r_v

    if (alpha < 0)
        fprintf('Alpha must be positive !');
        exit;
    end 
    gamma_x=1/dx^2;
    gamma_y=1/dy^2;
    
    % diagonal block matrix 
    A_diag=eye(N-1) *( 2*(gamma_x+gamma_y) + alpha +beta );
    A_diag =A_diag - diag(ones(N-2,1),1)*gamma_y ;
    A_diag =A_diag - diag(ones(N-2,1),-1) *gamma_y;

    
     % off diagonal block matrix 
    A_off= eye(N-1) *gamma_x ;   
    A=sparse((N-1)*(N-1), (N-1)*(N-1));
    
    % diagonal block matrices in A
    for i=1:N-1
    A((i-1)*(N-1)+1:(i-1)*(N-1)+(N-1), (i-1)*(N-1)+1:(i-1)*(N-1)+ (N-1))=A_diag;
    end 
   
    % off diagonal matrices in A
    for i=2:N-1
    A((i-2)*(N-1)+1:(i-2)*(N-1)+(N-1), (i-1)*(N-1)+1:(i-1)*(N-1)+ (N-1))=A_off; 
    A((i-1)*(N-1)+1:(i-1)*(N-1)+ (N-1),(i-2)*(N-1)+1:(i-2)*(N-1)+(N-1)) =A_off;
    end
    
    m =size (A);
    m =m(1);
    B=eye(m)*-beta;
    
    M=zeros(2*m);
    
    for ii = 1:2
        M((ii-1)*m+1:ii*m,(ii-1)*m+1:ii*m) = A;
        if ii~= 2
        M((ii-1)*m+1:ii*m,ii*m+1:(ii+1)*m) = B;
        M(ii*m+1:(ii+1)*m,(ii-1)*m+1:ii*m) = B;
        end
    end
    
    
