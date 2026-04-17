% Q5_Wave_FDM.m
clear; clc; close all;

% Setup: u_tt - c^2 u_xx = 0, t in [0, 0.5], x in [0, 0.5]
% c^2 = 1/(16*pi^2) -> c = 1/(4*pi)
c = 1/(4*pi); L = 0.5; T = 0.5;
f_init = @(x) 0;
g_init = @(x) sin(4*pi*x);
exact_sol = @(x, t) sin(t) .* sin(4*pi*x);

Nx = 20; h = L/Nx;
Nt = 40; k = T/Nt;
lambda = (c*k)/h; 
% Check CFL condition: lambda <= 1
if lambda > 1, warning('CFL condition violated!'); end

x = linspace(0, L, Nx+1)'; t = linspace(0, T, Nt+1);
U = zeros(Nx+1, Nt+1);
U(:,1) = f_init(x);

% First time step using initial velocity g(x)
for i = 2:Nx
    U(i,2) = (1 - lambda^2)*f_init(x(i)) + k*g_init(x(i)) ...
             + (lambda^2 / 2)*(f_init(x(i+1)) + f_init(x(i-1)));
end

% Explicit central difference scheme
for n = 2:Nt
    for i = 2:Nx
        U(i, n+1) = 2*(1 - lambda^2)*U(i,n) + lambda^2*(U(i+1,n) + U(i-1,n)) - U(i, n-1);
    end
end

[X, T_grid] = meshgrid(t, x);
U_exact = exact_sol(X, T_grid);

figure;
surf(T_grid, X, U); title('Q5(i): FDM Explicit Wave Equation');
xlabel('x'); ylabel('t'); zlabel('u(x,t)'); shading interp;