close all

N = 1.1;
T = 300;

Index = 1:10:T;


Tnum = length(Index);

MinEiga = zeros(Tnum,1);
MinEigc = zeros(Tnum,1);

Count = 1;
for k = 1:Tnum

    T = Index(k);
    A = zeros(T);
    C = zeros(T);
    for i = 1:T
        A(i,i) = T - i + 1;
        for j = i+1:T
            A(i,j) = sqrt((T-j+1)/N);
            C(i,j) = sqrt(1/N);
        end
    end

    A = A+A' - diag(diag(A));
    C = C+ C' + eye(T);

%     B = zeros(T);
%     B(1,:) = A(1,:);
%     B(:,1) = A(:,1);
    MinEiga(k) = min(eig(A));
    MinEigc(k) = min(eig(C));
end

histogram(MinEiga);
figure;histogram(MinEigc);

