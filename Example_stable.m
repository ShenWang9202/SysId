% ==========================================================================
% Unstable LTI systems for OLS estimators in the followingpaper
% Y. Zheng, N. Li, Non-asymptotic Identification of Linear Dynamical Systems  
%                                               using Multiple Trajectories
% ==========================================================================

clc;clear

% system data
n = 3;
m = 2;
p = 2;

A = randi(5,n,n); A = A./(1.2*max(abs(eig(A))));  % stable A
B = randi([-2,2],n,m);
C = randi([-2,2],p,n);
D = randi([-2,2],p,m);

T = 3*n;   % Length of Markov parameters % each experiments
G = D;
for k = 1:T-1
    G = [G, C*A^(k-1)*B];
end

% T1 = floor(T/2); T2 = T - T1 - 1;
% [hA,hB,hC,hD] = Ho_Kalman(G,T1,T2,n,m,p);

% noise level, and input singal
sigu = 1;
sigw = 0.1;
sigv = 0.1;


fprintf('Estimating an open-loop stable system ... \n\n')
fprintf('    Spectral radius of open-loop system:%6.3f \n', max(abs(eig(A))));

% --------------------------------------------------------------------
% Method 1: multiple trajectory
% --------------------------------------------------------------------

% Generating data 
N = 500;       % number of experiments
Y = zeros(p,N*T);
Z = zeros(m*T,N*T);
for i = 1:N    % generating trajectories
    [yi,ui,Zi]         = LTIsim(A,B,C,D,T,sigu,sigw,sigv);   
    Y(:,(i-1)*T+1:i*T) = yi;   % prepare the data for regression
    Z(:,(i-1)*T+1:i*T) = Zi;
end

hG1 = Y*pinv(Z);      % least squares solution

% -------------------------------------------------------------------- 
% Method 2: Single Trajectory
% --------------------------------------------------------------------
if max(abs(eig(A))) <= 1
    Nbar    = N*T;
    [y,u,~] = LTIsim(A,B,C,D,Nbar,sigu,sigw,sigv);  % generating data
    Y = y(:,T:end)';
    U = zeros(Nbar - T + 1,m*T);
    for i = 1:Nbar - T + 1
        ui = u(:,i+T-1:-1:i);
        U(i,:) = ui(:)';
    end
    %Uinv = (U'*U)^(-1)*U';
    hG2 = pinv(U)*Y;    % least-sqaure solutions
    hG2 = hG2';
end

% Estimation error
Gerror = [norm(hG1-G)./norm(G),norm(hG2-G)./norm(G)];  % spectral norm
fprintf('    The relative estimation errors of G: \n');
fprintf('            1) Mutiple trajectories: %6.3E\n',Gerror(1));
fprintf('            2) Single trajectory   : %6.3E\n',Gerror(2));

%% Ho-Kalman algorithm 
T1 = floor(T/2); T2 = T - T1 - 1;
flag = 0;       % a heuristic to choose system order or not
[hA1,hB1,hC1,hD1] = Ho_Kalman(hG1,T1,T2,n,m,p,flag);
[hA2,hB2,hC2,hD2] = Ho_Kalman(hG2,T1,T2,n,m,p);

% Hinf check 
ss0 = ss(A,B,C,D,[]);
ss2 = ss(hA1,hB1,hC1,hD1,[]);
ss3 = ss(hA2,hB2,hC2,hD2,[]);
Hinferror = [norm(ss2 - ss0,'inf'),norm(ss3 - ss0,'inf')]/norm(ss0,'inf'); 
fprintf('    The relative Hinf error after the Ho-Kalman algorithm: \n');
fprintf('            1) Mutiple trajectories: %6.3E\n',Hinferror(1));
fprintf('            2) Single trajectory   : %6.3E\n',Hinferror(2));

