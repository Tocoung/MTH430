% MTH 430 Lab Assignment 9 - Question 2: Nonlinear Shooting Method
clear; clc; close all;

% Given Parameters for BVP (i): y'' = -(y')^2 - y + ln(t), t in [1,2], y(1)=0, y(2)=ln(2)
a = 1; b = 2; alpha = 0; beta = log(2);

% f(t, y, y') = y''
f   = @(t, y, yp) -(yp.^2) - y + log(t);
fy  = @(t, y, yp) -1;
fyp = @(t, y, yp) -2*yp;

% Actual solution
exact_sol = @(t) log(t);

% --- Task (ii): Run with N = 2^3, Display Steps and Plot ---
N_plot = 2^3;
z_guess = 0.5; % Initial guess for y'(a)
tol = 1e-6;
max_iter = 20;

fprintf('--- Newton Iterations for Nonlinear Shooting (N = %d) ---\n', N_plot);
[t_num, y_num, iters] = nonlinear_shooting_newton(f, fy, fyp, a, b, alpha, beta, N_plot, z_guess, tol, max_iter, true);

% Fine mesh for actual solution
t_fine = linspace(a, b, 2^9);
y_fine = exact_sol(t_fine);

figure(3);
plot(t_fine, y_fine, 'b-', 'LineWidth', 1.5); hold on;
plot(t_num, y_num, 'ro--', 'LineWidth', 1.5);
title('Q2(i): Nonlinear Shooting Method vs Actual Solution');
xlabel('t'); ylabel('y(t)');
legend('Actual Solution (N=2^9)', 'Numerical Solution (N=2^3)', 'Location', 'best');
grid on; hold off;

% --- Task (iii): Error Analysis and Convergence ---
N_vals = [2^2, 2^3, 2^4, 2^5];
errors = zeros(size(N_vals));
h_vals = (b - a) ./ N_vals;

for i = 1:length(N_vals)
    [t_err, y_err, ~] = nonlinear_shooting_newton(f, fy, fyp, a, b, alpha, beta, N_vals(i), z_guess, tol, max_iter, false);
    y_exact_err = exact_sol(t_err);
    errors(i) = max(abs(y_err - y_exact_err));
end

figure(4);
loglog(h_vals, errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(h_vals, errors(1)*(h_vals/h_vals(1)).^4, 'r--', 'LineWidth', 1.5);
title('Q2: h vs Global Error (Log-Log Scale)');
xlabel('Step size (h)'); ylabel('Global Error (Max Norm)');
legend('Global Error', 'O(h^4) Reference', 'Location', 'best');
grid on; hold off;


% --- Helper Function: Nonlinear Shooting via Newton's Method ---
function [t, y, k] = nonlinear_shooting_newton(f, fy, fyp, a, b, alpha, beta, N, z0, tol, max_iter, display_steps)
    h = (b - a) / N;
    t = linspace(a, b, N+1);
    z_k = z0;
    
    if display_steps
        fprintf('Iter (k)\t z_k (Slope Guess)\t Error at y(b)\n');
        fprintf('---------------------------------------------------\n');
    end
    
    for k = 1:max_iter
        % System: 
        % u1 = y, u2 = y'
        % u3 = dy/dz, u4 = dy'/dz
        % F(t, U) = [u2; f(t,u1,u2); u4; fy*u3 + fyp*u4]
        U = zeros(4, N+1);
        U(:,1) = [alpha; z_k; 0; 1];
        
        system_eqs = @(t, u) [u(2); 
                              f(t, u(1), u(2)); 
                              u(4); 
                              fy(t, u(1), u(2))*u(3) + fyp(t, u(1), u(2))*u(4)];
                          
        % RK4 Integration
        for i = 1:N
            ti = t(i);
            k1 = system_eqs(ti, U(:,i));
            k2 = system_eqs(ti + h/2, U(:,i) + h/2 * k1);
            k3 = system_eqs(ti + h/2, U(:,i) + h/2 * k2);
            k4 = system_eqs(ti + h, U(:,i) + h * k3);
            U(:,i+1) = U(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        end
        
        y_b = U(1, end);
        z_b = U(3, end); % Derivative of y at b with respect to z
        
        error = y_b - beta;
        
        if display_steps
            fprintf('%d\t\t %.6f\t\t %.2e\n', k, z_k, abs(error));
        end
        
        if abs(error) < tol
            break;
        end
        
        % Newton update
        z_k = z_k - error / z_b;
    end
    
    y = U(1, :);
end