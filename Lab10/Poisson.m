% Q3_Poisson_FDM.m
clear; clc; close all;

% Setup: u_xx + u_yy = 0 on (0,1)x(0,1)
% BCs: u(x,0)=0, u(x,1)=x, u(0,y)=0, u(1,y)=y
L_x = 1; L_y = 1;
exact_sol = @(x, y) x .* y;

Nx = 32; Ny = 32;
hx = L_x/Nx; hy = L_y/Ny;
x = linspace(0, L_x, Nx+1); y = linspace(0, L_y, Ny+1);
[X, Y] = meshgrid(x, y);
U = zeros(Ny+1, Nx+1);

% Apply BCs
U(1, :) = 0;          % u(x,0) = 0
U(end, :) = x;        % u(x,1) = x
U(:, 1) = 0;          % u(0,y) = 0
U(:, end) = y';       % u(1,y) = y

% Inner nodes setup
Nx_in = Nx-1; Ny_in = Ny-1;
I_x = speye(Nx_in); I_y = speye(Ny_in);
D_x = spdiags(ones(Nx_in,1)*[1 -2 1], -1:1, Nx_in, Nx_in) / hx^2;
D_y = spdiags(ones(Ny_in,1)*[1 -2 1], -1:1, Ny_in, Ny_in) / hy^2;

% 2D Laplacian Matrix
A = kron(I_y, D_x) + kron(D_y, I_x);
F = zeros(Ny_in, Nx_in);

% Boundary condition forcing
F(1, :) = F(1, :) - U(1, 2:end-1) / hy^2;
F(end, :) = F(end, :) - U(end, 2:end-1) / hy^2;
F(:, 1) = F(:, 1) - U(2:end-1, 1) / hx^2;
F(:, end) = F(:, end) - U(2:end-1, end) / hx^2;

U_in = A \ F(:);
U(2:end-1, 2:end-1) = reshape(U_in, Ny_in, Nx_in);

figure;
surf(X, Y, U); title('Q3(i): FDM Solution for Poisson Equation');
xlabel('x'); ylabel('y'); zlabel('u(x,y)'); shading interp;