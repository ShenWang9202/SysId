% ==========================================================================
% An open-loop stable LTI system for OLS estimators in the following paper
% Y. Zheng, N. Li,Non-asymptotic Identification 
%                 of Linear Dynamical Systems using Multiple Trajectories
% ==========================================================================

clear

% system data
n = 3;
m = 2;
p = 2;

% unstable A; the first method is suitable for unstable systems as well
A = randi(5,n,n); A = A./(0.9*max(abs(eig(A))));  
B = randi([-2,2],n,m);
C = randi([-2,2],p,n);
D = randi([-2,2],p,m);

T = 3*n;   % Length of Markov parameters % each experiments
G = D;
for k = 1:T-1
    G = [G, C*A^(k-1)*B];
end

% noise level and input
sigu = 1;
sigw = 0.05;
sigv = 0.1;

fprintf('\nEstimating an open-loop unstable system ... \n\n')
fprintf('    Spectral radius of open-loop system:%6.3f \n', max(abs(eig(A))));

% --------------------------------------------------------------------
% Method 1: multiple trajectory 
% --------------------------------------------------------------------
N = 1000;  % number of experiments
Y = zeros(p,N*T);
Z = zeros(m*T,N*T);
for i = 1:N    % generating trajectories
    [yi,ui,Zi]         = LTIsim(A,B,C,D,T,sigu,sigw,sigv);   
    Y(:,(i-1)*T+1:i*T) = yi;   % prepare the data for regression
    Z(:,(i-1)*T+1:i*T) = Zi;
end

hG1 = Y*pinv(Z);      % least squares solution

% Estimation error
fprintf('    The relative estimation error of G:  %6.3E \n',norm(hG1-G)./norm(G));


%% Ho-Kalman algorithm 
T1 = floor(T/2); T2 = T - T1 - 1;
[hA1,hB1,hC1,hD1] = Ho_Kalman(hG1,T1,T2,n,m,p);

O  = obsv(A,C); 
O1 = obsv(hA1,hC1);
T  = (O'*O)^(-1)*O'*O1;     % find a similarity transformation
A1 = T*hA1*T^(-1); B1 = T*hB1; C1 = hC1*T^(-1);

Error = [norm(A-A1)/norm(A),norm(B-B1)/norm(B),norm(C-C1)/norm(C),norm(D-hD1)/norm(D)];
fprintf('    The relative estimation error of a state-space relaixation using Ho-Kalman algorithm \n');
fprintf('              |A - hA|/|A|: %6.3E \n',Error(1));
fprintf('              |B - hB|/|B|: %6.3E \n',Error(2));
fprintf('              |C - hC|/|C|: %6.3E \n',Error(3));
fprintf('              |D - hD|/|D|: %6.3E \n',Error(4));
