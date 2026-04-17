% Q2_NonlinearBVP_FDM.m
clear; clc; close all;

% Setup: y'' = -(y')^2 - y + ln(t), t in [1,2], y(1)=0, y(2)=ln(2)
a = 1; b = 2; alpha = 0; beta = log(2);
f = @(t, y, yp) -(yp.^2) - y + log(t);
fy = @(t, y, yp) -1; 
fyp = @(t, y, yp) -2*yp;
exact_sol = @(t) log(t);

N = 16; h = (b-a)/N; t = linspace(a, b, N+1)';
% Initial guess (linear interpolation)
y = alpha + (beta - alpha) * (t - a) / (b - a);
tol = 1e-6; max_iter = 20;

for iter = 1:max_iter
    F = zeros(N-1, 1); J = zeros(N-1, N-1);
    for i = 1:N-1
        ti = t(i+1);
        yi = y(i+1);
        yp = (y(i+2) - y(i)) / (2*h);
        
        F(i) = y(i+2) - 2*yi + y(i) - h^2 * f(ti, yi, yp);
        J(i,i) = -2 - h^2 * fy(ti, yi, yp);
        if i > 1
            J(i, i-1) = 1 + (h/2) * fyp(ti, yi, yp);
        end
        if i < N-1
            J(i, i+1) = 1 - (h/2) * fyp(ti, yi, yp);
        end
    end
    dy = J \ (-F);
    y(2:end-1) = y(2:end-1) + dy;
    if norm(dy, inf) < tol, break; end
end

figure;
plot(linspace(a,b,500), exact_sol(linspace(a,b,500)), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y, 'ro--', 'LineWidth', 1.5);
title('Q2(i): FDM for Nonlinear BVP'); xlabel('t'); ylabel('y(t)'); grid on;