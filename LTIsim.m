function [y,u,Z] = LTIsim(A,B,C,D,T,sigu,sigw,sigv)
    % Simulate the dynamical system 
    % x_{k+1} = Ax_k + Bu_k + Bw_k
    % y_k     = Cx_k + Du_k + v_k
    % u_k, w_k, v_k are Guassian vectors
    
    n = size(A,1);
    m = size(B,2);
    p = size(C,1);
    
    x = zeros(n,T);
    y = zeros(p,T);
    u = zeros(m,T);
    
    x(:,1) = zeros(n,1);
    u(:,1) = normrnd(0,sigu,[m,1]);
    w      = normrnd(0,sigw,[m,1]);
    v      = normrnd(0,sigv,[p,1]);
    x(:,2) = A*x(:,1) + B*u(:,1) + B*w;
    y(:,1) = C*x(:,1) + D*u(:,1) + v;    % initial output
    
    for k  = 2:T
        u(:,k) = normrnd(0,sigu,[m,1]);
        w        = normrnd(0,sigw,[m,1]);
        v        = normrnd(0,sigv,[p,1]);  
        
        if k < T
            x(:,k+1) = A*x(:,k) + B*u(:,k) + B*w;
        end
        
        y(:,k)   = C*x(:,k) + D*u(:,k) + v;   % output 
    end
   
    %% input matrix Zi
    Z = zeros(m*T,T);
    for k = 1:T
        Z((k-1)*m+1:k*m,k:end) = u(:,1:T-k+1);
    end
    
end