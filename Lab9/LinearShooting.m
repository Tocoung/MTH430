% MTH 430 Lab Assignment 9 - Question 1: Linear Shooting Method
clear; clc; close all;

% Given Parameters for BVP (i): y'' = 4(y-t), t in [0,1], y(0)=0, y(1)=2
% Standard form: y'' = p(t)y' + q(t)y + r(t)
a = 0; b = 1; alpha = 0; beta = 2;
p = @(t) 0;
q = @(t) 4;
r = @(t) -4*t;

% Actual solution for comparison
exact_sol = @(t) exp(2)/(exp(4)-1) * (exp(2*t) - exp(-2*t)) + t;

% --- Task (ii): Run with N = 2^3, Display Steps and Plot ---
N_plot = 2^3;
fprintf('--- Linear Shooting Method Steps for N = %d ---\n', N_plot);
[t_num, y_num] = linear_shooting_rk4(p, q, r, a, b, alpha, beta, N_plot, true);

% Fine mesh for actual solution
t_fine = linspace(a, b, 2^9);
y_fine = exact_sol(t_fine);

figure(1);
plot(t_fine, y_fine, 'b-', 'LineWidth', 1.5); hold on;
plot(t_num, y_num, 'ro--', 'LineWidth', 1.5);
title('Q1(i): Linear Shooting Method vs Actual Solution');
xlabel('t'); ylabel('y(t)');
legend('Actual Solution (N=2^9)', 'Numerical Solution (N=2^3)', 'Location', 'best');
grid on; hold off;

% --- Task (iii): Error Analysis and Convergence ---
N_vals = [2^2, 2^3, 2^4, 2^5];
errors = zeros(size(N_vals));
h_vals = (b - a) ./ N_vals;

for i = 1:length(N_vals)
    [t_err, y_err] = linear_shooting_rk4(p, q, r, a, b, alpha, beta, N_vals(i), false);
    y_exact_err = exact_sol(t_err);
    errors(i) = max(abs(y_err - y_exact_err));
end

figure(2);
loglog(h_vals, errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
% Reference slope for O(h^4) since we use RK4
loglog(h_vals, errors(1)*(h_vals/h_vals(1)).^4, 'r--', 'LineWidth', 1.5);
title('Q1: h vs Global Error (Log-Log Scale)');
xlabel('Step size (h)'); ylabel('Global Error (Max Norm)');
legend('Global Error', 'O(h^4) Reference', 'Location', 'best');
grid on; hold off;


% --- Helper Function: Linear Shooting with RK4 ---
function [t, y] = linear_shooting_rk4(p, q, r, a, b, alpha, beta, N, display_steps)
    h = (b - a) / N;
    t = linspace(a, b, N+1);
    
    % IVP 1: u'' = p*u' + q*u + r, u(a) = alpha, u'(a) = 0
    u = zeros(2, N+1); u(:,1) = [alpha; 0];
    % IVP 2: v'' = p*v' + q*v,     v(a) = 0,     v'(a) = 1
    v = zeros(2, N+1); v(:,1) = [0; 1];
    
    % RK4 System definitions
    f_u = @(t, u) [u(2); p(t)*u(2) + q(t)*u(1) + r(t)];
    f_v = @(t, v) [v(2); p(t)*v(2) + q(t)*v(1)];
    
    for i = 1:N
        ti = t(i);
        
        k1u = f_u(ti, u(:,i));
        k2u = f_u(ti + h/2, u(:,i) + h/2 * k1u);
        k3u = f_u(ti + h/2, u(:,i) + h/2 * k2u);
        k4u = f_u(ti + h, u(:,i) + h * k3u);
        u(:,i+1) = u(:,i) + h/6 * (k1u + 2*k2u + 2*k3u + k4u);
        
        k1v = f_v(ti, v(:,i));
        k2v = f_v(ti + h/2, v(:,i) + h/2 * k1v);
        k3v = f_v(ti + h/2, v(:,i) + h/2 * k2v);
        k4v = f_v(ti + h, v(:,i) + h * k3v);
        v(:,i+1) = v(:,i) + h/6 * (k1v + 2*k2v + 2*k3v + k4v);
    end
    
    % Combine solutions
    w2 = (beta - u(1, end)) / v(1, end);
    y = u(1, :) + w2 * v(1, :);
    
    if display_steps
        fprintf('Node (t)\t u_1(t)\t\t v_1(t)\t\t y(t)\n');
        for i = 1:N+1
            fprintf('%.4f\t\t %.6f\t %.6f\t %.6f\n', t(i), u(1,i), v(1,i), y(i));
        end
        fprintf('\n');
    end
end