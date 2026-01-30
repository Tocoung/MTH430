% MTH 430: Lab Assignment 4 - Consolidated Comparison
% Comparing CGS, MGS, Householder, and Givens for Least Squares
% -------------------------------------------------------------
clc; clear; close all;

% --- STEP 1: PREPARE DATA (Assignment 1, Table A) ---
t = [1; 3; 5; 6; 7; 8; 9; 10];
y = [1.3; 3.5; 4.2; 5.0; 7.0; 8.8; 10.1; 12.5];

% Define Polynomial Degree (e.g., 1 for Line, 2 for Quadratic)
degree = 7; 
N = length(t);

% Construct Design Matrix A (Vandermonde)
A_design = ones(N, degree + 1);
for p = 1:degree
    A_design(:, p+1) = t.^p;
end

fprintf('Analyzing Least Squares for Degree %d\n', degree);
[m_orig, n_orig] = size(A_design);

%% --- METHOD 1: CLASSICAL GRAM-SCHMIDT (Your Snippet 1) ---
fprintf('\n--- Method 1: Classical Gram-Schmidt ---\n');
A = A_design; % Reset A
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);

% Your Code:
for k = 1:n
    v = A(:,k);
    for j = 1:k-1
        R(j,k) = Q(:,j)' * A(:,k);
        v = v - R(j,k) * Q(:,j);
    end
    R(k,k) = norm(v);
    Q(:,k) = v / R(k,k);
end
Q_cgs = Q; R_cgs = R;

% Solve Rx = Q'y
b = Q_cgs' * y;
x_cgs = zeros(n, 1);
for i = n:-1:1
    sum_val = 0;
    for j = i+1:n
        sum_val = sum_val + R_cgs(i,j) * x_cgs(j);
    end
    x_cgs(i) = (b(i) - sum_val) / R_cgs(i,i);
end
disp('Coefficients (CGS):'); disp(x_cgs');


%% --- METHOD 2: MODIFIED GRAM-SCHMIDT (Your Snippet 2) ---
fprintf('\n--- Method 2: Modified Gram-Schmidt ---\n');
A = A_design; % Reset A
[m,n] = size(A);
Q = A;
R = zeros(n,n);

% Your Code:
for k = 1:n
    R(k,k) = norm(Q(:,k));
    Q(:,k) = Q(:,k) / R(k,k);
    for j = k+1:n
        R(k,j) = Q(:,k)' * Q(:,j);
        Q(:,j) = Q(:,j) - R(k,j) * Q(:,k);
    end
end
Q_mgs = Q; R_mgs = R;

% Solve Rx = Q'y
b = Q_mgs' * y;
x_mgs = zeros(n, 1);
for i = n:-1:1
    sum_val = 0;
    for j = i+1:n
        sum_val = sum_val + R_mgs(i,j) * x_mgs(j);
    end
    x_mgs(i) = (b(i) - sum_val) / R_mgs(i,i);
end
disp('Coefficients (MGS):'); disp(x_mgs');


%% --- METHOD 3: HOUSEHOLDER QR (Your Snippet 3) ---
fprintf('\n--- Method 3: Householder QR ---\n');
A = A_design; % Reset A
[m,n] = size(A);
R = A;
Q = eye(m);
Hs = {};

% Your Code (Wrapped to handle non-square A logic):
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
end
Q_hh = Q; R_hh = R;

% Solve Rx = Q'y
% Note: Householder produces Full QR (m x m Q). 
% We must truncate for Least Squares solution.
b_full = Q_hh' * y;
b = b_full(1:n);      % Keep top n elements
R_solve = R_hh(1:n, :); % Keep top n rows

x_hh = zeros(n, 1);
for i = n:-1:1
    sum_val = 0;
    for j = i+1:n
        sum_val = sum_val + R_solve(i,j) * x_hh(j);
    end
    x_hh(i) = (b(i) - sum_val) / R_solve(i,i);
end
disp('Coefficients (Householder):'); disp(x_hh');


%% --- METHOD 4: GIVENS QR (Previous Code) ---
fprintf('\n--- Method 4: Givens QR ---\n');
A = A_design; % Reset A
[m,n] = size(A);
R = A;
Q = eye(m);

% Givens Logic:
for k = 1:min(m-1,n)
    for j = k+1:m
        if R(j,k) ~= 0
            val_diag = R(k,k);
            val_elim = R(j,k);
            r = sqrt(val_diag^2 + val_elim^2);
            c = val_diag / r;
            s = -val_elim / r;
            
            G = eye(m);
            G(k,k) = c; G(j,j) = c;
            G(k,j) = -s; G(j,k) = s;
            
            R = G * R;
            Q = Q * G';
        end
    end
end
Q_giv = Q; R_giv = R;

% Solve Rx = Q'y
% Givens also produces Full QR (m x m Q)
b_full = Q_giv' * y;
b = b_full(1:n);      % Keep top n elements
R_solve = R_giv(1:n, :); % Keep top n rows

x_giv = zeros(n, 1);
for i = n:-1:1
    sum_val = 0;
    for j = i+1:n
        sum_val = sum_val + R_solve(i,j) * x_giv(j);
    end
    x_giv(i) = (b(i) - sum_val) / R_solve(i,i);
end
disp('Coefficients (Givens):'); disp(x_giv');


%% --- PLOTTING COMPARISON ---
t_plot = linspace(min(t), max(t), 100)';
y_cgs = zeros(size(t_plot));
y_mgs = zeros(size(t_plot));
y_hh  = zeros(size(t_plot));
y_giv = zeros(size(t_plot));

for p = 0:degree
    y_cgs = y_cgs + x_cgs(p+1) * t_plot.^p;
    y_mgs = y_mgs + x_mgs(p+1) * t_plot.^p;
    y_hh  = y_hh  + x_hh(p+1)  * t_plot.^p;
    y_giv = y_giv + x_giv(p+1) * t_plot.^p;
end

figure;
plot(t, y, 'ko', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Data'); hold on;
plot(t_plot, y_cgs, 'r-', 'LineWidth', 2, 'DisplayName', 'Classical GS');
plot(t_plot, y_mgs, 'g--', 'LineWidth', 2, 'DisplayName', 'Modified GS');
plot(t_plot, y_hh,  'b:', 'LineWidth', 2, 'DisplayName', 'Householder');
plot(t_plot, y_giv, 'm-.', 'LineWidth', 2, 'DisplayName', 'Givens');

title('Comparison of QR Methods for Least Squares');
xlabel('t'); ylabel('y');
legend('Location', 'best');
grid on;
hold off;