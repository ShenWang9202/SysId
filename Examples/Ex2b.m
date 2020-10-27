% ==========================================================================
% Code for Example IV.B in our paper
% Y. Zheng, N. Li, Non-asymptotic Identification of Linear Dynamical Systems  
%                                               using Multiple Trajectories
% ==========================================================================

clc;clear
%% system data
n = 3;p = 1;m = 3;
A = [1 0.01 0; 0.01 1 0.01; 0 0.01 1];
B = eye(n);
C = [1,0,0];
D = zeros(p,m);


%% noise 
sigu = 1;
sigv = 0.5;

%% Generating data 
Num = 500;  % number of experiments

SigW = 0:0.2:0.8;
Tind = 10:5:40;

NumRepeat = 20;
Err1 = cell(length(SigW),1);
Err2 = cell(length(SigW),1);
Err3 = cell(length(SigW),1);

for indsig = 1: length(SigW)

    Err1{indsig} = zeros(NumRepeat,length(Tind));
    Err2{indsig} = zeros(NumRepeat,length(Tind));
    Err3{indsig} = zeros(NumRepeat,length(Tind));

    sigw = SigW(indsig);
    for Idx = 1:length(Tind)
        fprintf('Length of each experiment: %d\n',Tind(Idx));
        for Re = 1:NumRepeat
            fprintf('   Number of trials: %d\n',Re);

            N = Num;
            T = Tind(Idx);
            Y = zeros(p,N*T);
            Z = zeros(m*T,N*T);

            G = D;
            for k = 1:T-1
                G = [G, C*A^(k-1)*B];
            end

            Y1 = zeros(p,N);
            Z1 = zeros(m*T,N);

            for i = 1:N
                [yi,ui,Zi] = LTIsim(A,B,C,D,T,sigu,sigw,sigv);

                Y(:,(i-1)*T+1:i*T) = yi;
                Z(:,(i-1)*T+1:i*T) = Zi;

                Y1(:,i) = yi(:,end);
                Z1(:,i) = Zi(:,end);

            end

       %% Method 1: Multiple Trajectory -- all data point
            hG1 = Y*pinv(Z);

       %% Method 2: Multiple Trajectory -- final data point
        hG2 = Y1*pinv(Z1);
               
       %% Method 3: Single Trajectory + prefilter
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
           Err1{indsig}(Re,Idx) = [norm(hG1-G)]./norm(G);
           Err2{indsig}(Re,Idx) = [norm(hG2-G)]./norm(G);
           Err3{indsig}(Re,Idx) = [norm(hG3-G)]./norm(G);
        end
    end
    save data_unstable_varyingT
end

%% Figure
paperFig4