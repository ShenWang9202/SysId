

T = 3;
n = 1;

A = randi([-5,5],n*T,n*T);
A = (A + A')/2;
%A = A - 1.1*min(eig(A))*eye(n*T);

B = zeros(T);
for i = 1:T
    B(i,i) = min(eig(A((i-1)*n+1:i*n,(i-1)*n+1:i*n)));
    for j = i + 1:T
        B(i,j) = max(svd(A((i-1)*n+1:i*n,(j-1)*n+1:j*n)));
    end
end

B = B + B' - diag(diag(B));

[min(eig(A)), min(eig(B))]