clc;clear
%% system data
n = 3;
m = 2;
p = 2;

A = randi(5,n,n); A = A./(0.5*max(abs(eig(A))));  % stable A
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
sigw = 0.1;
sigv = 0.1;

%% Method 1: multiple trajectory
% Generating data 
N = 50000;  % number of experiments

Y = zeros(p,N*T);
Z = zeros(m*T,N*T);

for i = 1:N
    [yi,ui,Zi] = LTIsim(A,B,C,D,T,sigu,sigw,sigv);
    
    Y(:,(i-1)*T+1:i*T) = yi;
    Z(:,(i-1)*T+1:i*T) = Zi;
    
end

% least squares solution
hG1 = Y*pinv(Z);

%% Method 2: Single Trajectory
if max(abs(eig(A))) <= 1
    Nbar    = N*T;
    [y,u,~] = LTIsim(A,B,C,D,Nbar,sigu,sigw,sigv);

    Y = y(:,T:end)';
    U = zeros(Nbar - T + 1,m*T);
    for i = 1:Nbar - T + 1
        ui = u(:,i+T-1:-1:i);
        U(i,:) = ui(:)';
    end
    Uinv = (U'*U)^(-1)*U';
    hG2 = Uinv*Y;
    hG2 = hG2';
end
%% Estimation error
norm(hG1-G)./norm(G)



%% Ho-Kalman algorithm 
T1 = floor(T/2); T2 = T - T1 - 1;
[hA1,hB1,hC1,hD1] = Ho_Kalman(hG1,T1,T2,n,m,p);
%[hA2,hB2,hC2,hD2] = Ho_Kalman(hG2,T1,T2,n,m,p);

O  = obsv(A,C); 
O1 = obsv(hA1,hC1);

T = (O'*O)^(-1)*O'*O1;
A1 = T*hA1*T^(-1); B1 = T*hB1; C1 = hC1*T^(-1);


[norm(A-A1)/norm(A),norm(B-B1)/norm(B),norm(C-C1)/norm(C)]

% O2 = obsv(hA2,hC2);
% T = (O'*O)^(-1)*O'*O1;
% A2 = T*hA2*T^(-1); B2 = T*hB2; C2 = hC2*T^(-1);
% 
% [B,B1,B2]
% [C;C1;C2]
% 
% norm(A-A2)/norm(A),norm(B-B2)/norm(B),norm(C-C2)/norm(C)]


% Hinf check 
% ss0 = ss(A,B,C,D,[]);
% ss1 = ss(hA,hB,hC,hD,[]);
% ss2 = ss(hA1,hB1,hC1,hD1,[]);
% ss3 = ss(hA2,hB2,hC2,hD2,[]);
% 
% [norm(ss1 - ss0,'inf'),norm(ss2 - ss0,'inf'),norm(ss3 - ss0,'inf')]/norm(ss0,'inf')
% 
% 
