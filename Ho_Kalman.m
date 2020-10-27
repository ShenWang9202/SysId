function [A,B,C,D,sigma] = Ho_Kalman(G,T1,T2,n,m,p,flag)
% A naive implementation of the Ho-Kalman algorithm to find (A,B,C,D)
% Input: 
%     G: markov paramter
% T1,T2: dimension of Hankel matrix  -- T = T1 + T2 + 1
% n,m,p: system dimension
% Output: State space realization A, B, C, D
    
    if nargin <=6
        flag = 0;  % a naive adaptation to choose system order.
    end
    
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
    if flag == 0
        for i = 1:n
            hS(i,i) = S(i,i);
        end
        L = U*hS*V';
    else
        totalSig = sum(diag(S));
        tmp = 0;
        for i = 1:min(size(S))
            hS(i,i) = S(i,i);
            tmp = tmp + hS(i,i);
            if tmp/totalSig > 0.9
                n = i;   % system order
                break;
            end
        end
        L = U*hS*V';
    end
      
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
    A = pinv(hO)*Hplu*pinv(hQ);
    %hOi  = (hO'*hO)^(-1)*hO';
    %hQi  = hQ'*(hQ*hQ')^(-1);   % psedoinverse
    %A = hOi*Hplu*hQi;
    
    D = G(:,1:m);  
end