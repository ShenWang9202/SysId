clc;clear
%% system data
n = 3;
m = 2;
p = 2;

% unstable A; the first method is suitable for unstable systems as well
A = randi(5,n,n); A = A./(0.8*max(abs(eig(A))));  
B = randi([-2,2],n,m);
C = randi([-2,2],p,n);
D = randi([-2,2],p,m);

T = 10;   % Length of Markov parameters % each experiments
G = D;
for k = 1:T-1
    G = [G, C*A^(k-1)*B];
end


sigu = 1;
sigw = 0.1;
sigv = 0.1;

%% Method 1: multiple trajectory
% Generating data 
N = 5000;  % number of experiments

Y = zeros(p,N*T);
Z = zeros(m*T,N*T);

for i = 1:N
    [yi,ui,Zi] = LTIsim(A,B,C,D,T,sigu,sigw,sigv);
    
    Y(:,(i-1)*T+1:i*T) = yi;
    Z(:,(i-1)*T+1:i*T) = Zi;
    
end

% least squares solution
hG1 = Y*pinv(Z);


%% Estimation error
norm(hG1-G)./norm(G)



%% Ho-Kalman algorithm 
T1 = floor(T/2); T2 = T - T1 - 1;
[hA1,hB1,hC1,hD1] = Ho_Kalman(hG1,T1,T2,n,m,p);

O  = obsv(A,C); 
O1 = obsv(hA1,hC1);

T = (O'*O)^(-1)*O'*O1;
A1 = T*hA1*T^(-1); B1 = T*hB1; C1 = hC1*T^(-1);


[norm(A-A1)/norm(A),norm(B-B1)/norm(B),norm(C-C1)/norm(C)]
