function [A,B,C,D] = Ho_Kalman(G,T,T1,T2,n,m,p)
    % Ho-Kalman algorithm to find (A,B,C,D)
    % T = T1 + T2 + 1
    
    % not work, June 08, 2020
    
    % step 1: Hankel matrix
    H = zeros(p*T1,m*(T2+1));  % Hankel matrix
    for i = 1:T1
        for j = 1:T2+1
            H((i-1)*p+1:i*p,(j-1)*m+1:j*m) = G(:,(i+j-1)*m+1:(i+j)*m);
        end
    end
    
    % step 2: first mT2 columns of H
    Hneg = H(:,1:m*T2);
    
    % step 3: rank-n-approximation of hH
    [U,S,V] = svd(Hneg);
    hS = zeros(size(S));
    for i = 1:n
        hS(i,i) = S(i,i);
    end
    L = U*hS*V';
      
    % SVD decomposition
    [U,Sig,V] = svd(L);
    Sigma = Sig(1:n,1:n);  %% only keep the non-zero elements
    Uc = U(:,1:n);
    Vc = V(:,1:n);
    
    hO = Uc*Sigma^(0.5);
    hQ = Sigma^(0.5)*Vc';
    
    C = hO(1:p,:);
    B = hQ(:,1:m);
    
    Hplu = H(:,m+1:end);
    hOi  = (hO'*hO)^(-1)*hO';
    hQi  = hQ'*(hQ*hQ')^(-1);   % psedoinverse
    A = hOi*Hplu*hQi;
    
    D = G(:,1:m);
     
    
end