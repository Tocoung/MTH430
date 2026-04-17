% Q4_Heat_FDM.m
clear; clc; close all;

% Setup: u_t = u_xx, t in [0, 0.1], x in [0, 2]
% u(t,0)=0, u(t,2)=0, u(0,x)=sin(pi*x/2)
L = 2; T = 0.1;
exact_sol = @(x, t) exp(-(pi^2/4)*t) .* sin(pi*x/2);

Nx = 20; Nt = 100; % Ensure stable k/h^2 ratio for explicit
h = L/Nx; k = T/Nt;
x = linspace(0, L, Nx+1)'; t = linspace(0, T, Nt+1);
lambda = k / h^2;

U_CN = zeros(Nx+1, Nt+1);
U_CN(:,1) = sin(pi*x/2);

% Crank-Nicolson Assembly (Theta = 0.5)
theta = 0.5;
A = spdiags([-lambda*theta*ones(Nx-1,1), (1+2*lambda*theta)*ones(Nx-1,1), -lambda*theta*ones(Nx-1,1)], -1:1, Nx-1, Nx-1);
B = spdiags([lambda*(1-theta)*ones(Nx-1,1), (1-2*lambda*(1-theta))*ones(Nx-1,1), lambda*(1-theta)*ones(Nx-1,1)], -1:1, Nx-1, Nx-1);

for n = 1:Nt
    rhs = B * U_CN(2:end-1, n);
    U_CN(2:end-1, n+1) = A \ rhs;
end

[X, T_grid] = meshgrid(t, x);
U_exact = exact_sol(X, T_grid);

figure;
subplot(1,2,1); surf(X, T_grid, U_CN); title('Crank-Nicolson'); xlabel('t'); ylabel('x'); shading interp;
subplot(1,2,2); surf(X, T_grid, U_exact); title('Actual Solution'); xlabel('t'); ylabel('x'); shading interp;