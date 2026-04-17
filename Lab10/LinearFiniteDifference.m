% Q1_LinearBVP_FDM.m
clear; clc; close all;

% Problem Setup: y'' = 4(y-t), t in [0,1], y(0)=0, y(1)=2
% Standard form: y'' = p(t)y' + q(t)y + r(t)
a = 0; b = 1; alpha = 0; beta = 2;
p = @(t) 0; q = @(t) 4; r = @(t) -4*t;
exact_sol = @(t) exp(2)/(exp(4)-1) * (exp(2*t) - exp(-2*t)) + t;

% --- (b) Numerical vs Actual Solution ---
N = 32; % Adjust N as needed
h = (b-a)/N; t = linspace(a, b, N+1)';
A = zeros(N-1, N-1); F = zeros(N-1, 1);

for i = 1:N-1
    ti = t(i+1);
    A(i,i) = -(2 + h^2 * q(ti));
    if i > 1, A(i, i-1) = 1 + (h/2)*p(ti); end
    if i < N-1, A(i, i+1) = 1 - (h/2)*p(ti); end
    F(i) = h^2 * r(ti);
end
F(1) = F(1) - (1 + (h/2)*p(t(2))) * alpha;
F(end) = F(end) - (1 - (h/2)*p(t(end-1))) * beta;

y_num = [alpha; A\F; beta];
t_fine = linspace(a, b, 500);

figure;
plot(t_fine, exact_sol(t_fine), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_num, 'ro--', 'LineWidth', 1.5);
title('Q1(i): FDM for Linear BVP'); xlabel('t'); ylabel('y(t)');
legend('Actual Solution', sprintf('Numerical (N=%d)', N), 'Location', 'best'); grid on;

% --- (c) Global Error vs h ---
N_vals = [8, 16, 32, 64, 128];
errors = zeros(size(N_vals));
h_vals = (b-a)./N_vals;

for k = 1:length(N_vals)
    Nk = N_vals(k); hk = h_vals(k); tk = linspace(a, b, Nk+1)';
    Ak = zeros(Nk-1, Nk-1); Fk = zeros(Nk-1, 1);
    for i = 1:Nk-1
        ti = tk(i+1);
        Ak(i,i) = -(2 + hk^2 * q(ti));
        if i > 1, Ak(i, i-1) = 1 + (hk/2)*p(ti); end
        if i < Nk-1, Ak(i, i+1) = 1 - (hk/2)*p(ti); end
        Fk(i) = hk^2 * r(ti);
    end
    Fk(1) = Fk(1) - (1 + (hk/2)*p(tk(2))) * alpha;
    Fk(end) = Fk(end) - (1 - (hk/2)*p(tk(end-1))) * beta;
    yk = [alpha; Ak\Fk; beta];
    errors(k) = max(abs(yk - exact_sol(tk)));
end

figure;
loglog(h_vals, errors, 'bo-', 'LineWidth', 2); hold on;
loglog(h_vals, errors(1)*(h_vals/h_vals(1)).^2, 'r--', 'LineWidth', 1.5);
title('Convergence: Global Error vs h'); xlabel('h'); ylabel('Max Error');
legend('FDM Error', 'O(h^2) Reference', 'Location', 'best'); grid on;