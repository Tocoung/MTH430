A=[1,0,0;
    2,-2,3;
    -1,1,0;
    0,0,1;
    1,0,-1];
[m,n] = size(A);
R = A;
Q = eye(m);
Hs = {};

for k = 1:min(m-1,n)
    x = R(k:m,k);
    e1 = zeros(length(x),1); e1(1) = 1;
    alpha = -sign(x(1)) * norm(x);
    if alpha == 0
        v = x;
    else
        v = x - alpha*e1;
    end
    if norm(v) == 0
        Hk_sub = eye(length(x));
    else
        v = v / norm(v);
        Hk_sub = eye(length(x)) - 2*(v*v');
    end
    Hk = eye(m);
    Hk(k:m,k:m) = Hk_sub;
    R = Hk * R;
    Q = Q * Hk';
    Hs{end+1} = Hk;
    disp(Hk);
end
disp(R);
disp(Q);