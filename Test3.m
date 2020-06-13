
m = 3;
T = 5;
Nt = 20:10:50

Repeat = 10;

minEig = zeros(Repeat,length(Nt));

for Index = 1:length(Nt)
    N = Nt(Index);
    for Indexj = 1:Repeat
        Z = zeros(m*T,N*T);
        W = zeros(m*T,N*T);
        for i = 1:N
            u = zeros(m,T);
            w = zeros(m,T);
            for k  = 1:T
                u(:,k) = normrnd(0,sigu,[m,1]);
                w(:,k) = normrnd(0,sigu,[m,1]);
            end 
            %% input matrix Zi
            Zi = zeros(m*T,T);
            Wi = zeros(m*T,T);
            for k = 1:T
                Zi((k-1)*m+1:k*m,k:end) = u(:,1:T-k+1);
                Wi((k-1)*m+1:k*m,k:end) = w(:,1:T-k+1);
            end

             Z(:,(i-1)*T+1:i*T) = Zi;
             W(:,(i-1)*T+1:i*T) = Wi;
        end
        minEig(Indexj,Index) = norm(W*Z',2);
    end
end

4*sqrt(Nt*2*m*T)

minEig

% 
% sqrt(Nt) - sqrt(m)
% (sqrt(Nt) - sqrt(m))*sqrt(T)