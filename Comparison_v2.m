clc;clear
%% system data
n = 4;
m = 2;
p = 2;

A = randi(5,n,n); A = A./(1.2*max(abs(eig(A))));  % stable A
B = randi([-2,2],n,m);
C = randi([-2,2],p,n);
D = randi([-2,2],p,m);

T = 10;   % Length of Markov parameters % each experiments
G = D;
for k = 1:T-1
    G = [G, C*A^(k-1)*B];
end

%  T1 = floor(T/2); T2 = T - T1 - 1;
%  [hA,hB,hC,hD] = Ho_Kalman(G,T1,T2,n,m,p);


sigu = 1;
sigw = 0.1; %0.05;
sigv = 0.2;

%% Method 1: multiple trajectory -- all outputs
% Generating data 

Num = 100:50:500;  % number of experiments

NumRepeat =10;

Err1 = zeros(NumRepeat,length(Num));
Err2 = zeros(NumRepeat,length(Num));
Err3 = zeros(NumRepeat,length(Num));

for Idx = 1:length(Num)
    fprintf('Number of samples: %d\n',Num(Idx));
    for Re = 1:NumRepeat
        fprintf('   Number of trials: %d\n',Re);
        
        N = Num(Idx);
        Y = zeros(p,N*T);
        Z = zeros(m*T,N*T);

        Y1 = zeros(p,N);
        Z1 = zeros(m*T,N);

        for i = 1:N
            [yi,ui,Zi] = LTIsim(A,B,C,D,T,sigu,sigw,sigv);

            Y(:,(i-1)*T+1:i*T) = yi;
            Z(:,(i-1)*T+1:i*T) = Zi;

            Y1(:,i) = yi(:,end);
            Z1(:,i) = Zi(:,end);

        end

        % least squares solution
        hG1 = Y*pinv(Z);

        %% Method 2: Multiple Trajectory -- final data point
        hG2 = Y1*pinv(Z1);


        %% Method 3: Single Trajectory
        if max(abs(eig(A))) <= 1
            Nbar    = N*T;
            [y,u,~] = LTIsim(A,B,C,D,Nbar,sigu,sigw,sigv);

            Y = y(:,T:end)';
            U = zeros(Nbar - T + 1,m*T);
            for i = 1:Nbar - T + 1
                ui = u(:,i+T-1:-1:i);
                U(i,:) = ui(:)';
            end
            %Uinv = (U'*U)^(-1)*U';
            hG3 = pinv(U)*Y;
            hG3 = hG3';
        else
            hG3 = G;
        end
        %% Estimation error
       Err1(Re,Idx) = [norm(hG1-G)]./norm(G);
       Err2(Re,Idx) = [norm(hG2-G)]./norm(G);
       Err3(Re,Idx) = [norm(hG3-G)]./norm(G);
    end
end

figure
shadedErrorBar(Num*T, Err1, {@mean,@std}, 'lineprops', '-r');
hold on
shadedErrorBar(Num*T, Err2, {@mean,@std}, 'lineprops', '-b')
hold on
shadedErrorBar(Num*T, Err3, {@mean,@std}, 'lineprops', '-m')