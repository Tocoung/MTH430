%% Gram–Schmidt (modified) — beginner friendly demo
% Columns of V are the input vectors v1, v2, ...
% Outputs:
%   U : orthogonal vectors (columns)
%   Q : orthonormal vectors (columns)

clc; clear; close all;

% -----------------------------
% Example input (from your problem)
% -----------------------------
% v1 = [2; -1], v2 = [1; 3]
V = [2, 1;
     -1, 3];

% -----------------------------
% Prepare storage
% -----------------------------
[m, n] = size(V);    % m = vector length, n = number of vectors
U = zeros(m, n);     % will store orthogonal vectors
Q = zeros(m, n);     % will store orthonormal vectors

% -----------------------------
% Modified Gram–Schmidt
% -----------------------------
for j = 1:n
    v = V(:, j);               % current input vector

    % subtract projections onto previously computed orthonormal vectors
    for i = 1:(j-1)
        r = Q(:, i)' * v;      % projection coefficient (scalar)
        v = v - r * Q(:, i);   % remove component along Q(:,i)
    end

    U(:, j) = v;               % orthogonal vector (possibly zero if dependent)

    % normalize to get orthonormal vector (if nonzero)
    nv = norm(v);
    if nv < 1e-12
        warning('Vector %d is (nearly) linearly dependent and becomes zero.', j);
        Q(:, j) = zeros(m,1);
    else
        Q(:, j) = v / nv;
    end
end

% -----------------------------
% Display results
% -----------------------------
disp('Input matrix V (columns are input vectors):');
disp(V);

disp('Orthogonal vectors U (columns):');
disp(U);

disp('Orthonormal vectors Q (columns):');
disp(Q);

% Check orthogonality / orthonormality (small numbers ~ 0 expected)
disp('Q''*Q  (should be identity for independent vectors):');
disp(Q' * Q);

% End of script
