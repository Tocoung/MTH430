A=[1,2,3;
    4,6,0;
    5,1,1;
    0,0,1];
[m,n] = size(A);
Q = A;
R = zeros(n,n);

for k = 1:n
    R(k,k) = norm(Q(:,k));
    Q(:,k) = Q(:,k) / R(k,k);
    for j = k+1:n
        R(k,j) = Q(:,k)' * Q(:,j);
        Q(:,j) = Q(:,j) - R(k,j) * Q(:,k);
    end
end
disp(Q);
disp(R);