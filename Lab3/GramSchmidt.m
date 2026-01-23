A=[0,20,-14;
    3,27,-4;
    4,11,-2];
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);

for k = 1:n
    v = A(:,k);
    for j = 1:k-1
        R(j,k) = Q(:,j)' * A(:,k);
        v = v - R(j,k) * Q(:,j);
    end
    R(k,k) = norm(v);
    Q(:,k) = v / R(k,k);
end
disp(Q);
disp(R);