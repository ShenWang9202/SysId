clc;clear
%% system data
n = 4;
m = 2;
p = 1;

A = randi(10,n,n); A = A./(0.5*max(abs(eig(A))));  % stable A
B = randi([-2,2],n,m);
C = randi([-2,2],p,n);
D = randi([-2,2],p,m);

%% Generating data 
N = 2000;  % number of experiments
T = 10;   % Length of each experiment

G = D;
for k = 1:T-1
    G = [G, C*A^(k-1)*B];
end

Y = zeros(p,N*T);
Z = zeros(m*T,N*T);

sigu = 1;
sigw = 0.1;
sigv = 0.1;
for i = 1:N
    [yi,ui,Zi] = LTIsim(A,B,C,D,T,sigu,sigw,sigv);
    
    Y(:,(i-1)*T+1:i*T) = yi;
    Z(:,(i-1)*T+1:i*T) = Zi;
    
end


%% least squares solution
hG = Y*pinv(Z);

[G;
hG;
hG - G]
norm(hG - G)%./norm(G)
max(abs(eig(A)))
%% Ho-Kalman algorithm 
% T1 = floor(T/2); T2 = T - T1 - 1;
% [hA,hB,hC,hD] = Ho_Kalman(G,T,T1,T2,n,m,p);
% 
% %% balanced realization
% ss1 = ss(A,B,C,D,[]);
% [ssb,~] = balreal(ss1);
% 
% hA - ssb.A
% hB - ssb.B
% hC - ssb.C
% hD - ssb.D
% 

