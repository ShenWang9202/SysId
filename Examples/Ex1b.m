% =================================================================
%
% Code for Example IV.A in our paper
% Y. Zheng, N. Li, Non-asymptotic  Identification  of  Partially  Observable  
%                     Linear Time-invariant  Systems  using  Multiple  Trajectories
%
% =================================================================

clc;clear

%% Example 2: unstable system
n = 3;p = 1;m = 3;
A = [1 0.01 0; 0.01 1 0.01; 0 0.01 1];
B = eye(n);
C = [1,0,0];
D = zeros(p,m);

T = 10;   % Length of Markov parameters

G = D;
for k = 1:T-1
    G = [G, C*A^(k-1)*B];
end

sigu = 1;    % noise levels
sigw = 0.2;
sigv = 0.5;

% Generating data 
Num = 50:50:500;  % number of experiments

%% Method 1: multiple trajectory -- all outputs

NumRepeat =20;       % each experiment repeat 20 times
Err1 = zeros(NumRepeat,length(Num));
Err2 = zeros(NumRepeat,length(Num));
Err3 = zeros(NumRepeat,length(Num));

for Idx = 1:length(Num)
    fprintf('Number of samples: %d\n',Num(Idx));
     N = Num(Idx);
    for Re = 1:NumRepeat
        fprintf('   Number of trials: %d\n',Re);
        Y = zeros(p,N*T);
        Z = zeros(m*T,N*T);

        Y1 = zeros(p,N);
        Z1 = zeros(m*T,N);
        for i = 1:N
            [yi,ui,Zi] = LTIsim(A,B,C,D,T,sigu,sigw,sigv);

            Y(:,(i-1)*T+1:i*T) = yi;
            Z(:,(i-1)*T+1:i*T) = Zi;

            Y1(:,i) = yi(:,end);  % final data point
            Z1(:,i) = Zi(:,end);
        end

        % ====== Method 1: Multiple Trajectory -- all data point =======
        hG1 = Y*pinv(Z);

        % ====== Method 2: Multiple Trajectory -- final data point =====
        hG2 = Y1*pinv(Z1);

        % ====== Method 3: Single Trajectory  + prefilter   ============
        if max(abs(eig(A))) <= 1
            Nbar    = N*T;
            [y,u,~] = LTIsim(A,B,C,D,Nbar,sigu,sigw,sigv);
            
            L = T*1;
            mu = 10;
            N1 = L*T+1;
            K = zeros(p*L,Nbar - N1 + 1);
            Y = y(:,N1:end);
            U = zeros(m*T,Nbar - N1 + 1);
            for i = N1:Nbar
                ui = u(:,i:-1:i-T+1);
                U(:,i-N1+1) = ui(:);
                kt = y(:,i-T:-T:i-T*L);
                K(:,i-N1+1) = kt(:);
            end
            Phi = Y*(K'*(K*K'+2*mu^2*eye(L*p))^(-1));
            hG3 = (Y-Phi*K)*pinv(U);       
        else
            hG3 = G;
        end
        %% Estimation error
       Err1(Re,Idx) = [norm(hG1-G)]./norm(G);
       Err2(Re,Idx) = [norm(hG2-G)]./norm(G);
      % Err3(Re,Idx) = [norm(hG3-G)]./norm(G);  %this method do not work
    end
end

save data_unstable
% draw figure
paperFig2