clc;
clear;
close all;

% ---------------------------------
% Given observational data
% ---------------------------------
x = [1.02; 0.95; 0.87; 0.77; 0.67; 0.56; 0.44; 0.30; 0.16; 0.01];
y = [0.39; 0.32; 0.27; 0.22; 0.18; 0.15; 0.13; 0.12; 0.13; 0.15];

% ---------------------------------
% Construct the design matrix A
% Model: a*y^2 + b*x*y + c*x + d*y + e = x^2
% ---------------------------------
A = [y.^2, x.*y, x, y, ones(size(x))];

% Right-hand side vector
b = x.^2;

% ---------------------------------
% Normal equations
% ---------------------------------
ATA = A' * A;
ATb = A' * b;

% ---------------------------------
% Cholesky decomposition
% ATA = L * L'
% ---------------------------------
L = chol(ATA, 'lower');

% Forward substitution: L*z = ATb
z = L \ ATb;

% Back substitution: L'*p = z
p = L' \ z;

% Extract parameters
a = p(1);
b_coef = p(2);
c = p(3);
d = p(4);
e = p(5);

% Display results
fprintf('Least squares parameters:\n');
fprintf('a = %.6f\n', a);
fprintf('b = %.6f\n', b_coef);
fprintf('c = %.6f\n', c);
fprintf('d = %.6f\n', d);
fprintf('e = %.6f\n', e);

% ---------------------------------
% Plot the orbit and data
% ---------------------------------
y_plot = linspace(min(y)-0.05, max(y)+0.05, 500);

% Quadratic coefficients in x
B = b_coef*y_plot + c;
C = a*y_plot.^2 + d*y_plot + e;

% Discriminant
D = B.^2 + 4*C;

% Only real solutions
valid = D >= 0;

x1 = (B(valid) + sqrt(D(valid))) / 2;
x2 = (B(valid) - sqrt(D(valid))) / 2;
y_valid = y_plot(valid);

figure;
hold on;
grid on;

% Plot observed data
plot(x, y, 'ro', 'MarkerSize', 7, 'LineWidth', 2);

% Plot fitted orbit (two branches)
plot(x1, y_valid, 'b', 'LineWidth', 2);
plot(x2, y_valid, 'b', 'LineWidth', 2);

xlabel('x');
ylabel('y');
title('Least Squares Elliptical Orbit Fit');
legend('Observed data', 'Fitted orbit');
axis equal;
