function euler_method_ivp(caseId)
%EULER_METHOD_IVP  Forward Euler solver for 4 given IVPs.
%   euler_method_ivp(caseId)
%   caseId = 1,2,3,4 selects the ODE from the assignment.
%
%   This script:
%   1) prints iteration values for N = 2^3,
%   2) plots exact solution (N = 2^9) vs Euler solution (N = 2^3),
%   3) plots global error vs h for N = 2, 2^2, 2^3, 2^4.

if nargin < 1
    caseId = 1;
end

prob = getProblem(caseId);

% 1) Show iteration steps on a coarse grid
Ncoarse = 2^3;
[tc, yc] = eulerSolve(prob.f, prob.t0, prob.tf, prob.y0, Ncoarse);
printIterations(tc, yc, prob.exact, sprintf('Forward Euler | Case %d | N = %d', caseId, Ncoarse));

% 2) Compare coarse numerical solution with a fine exact mesh
Nfine = 2^9;
[tf, yexactFine] = exactOnGrid(prob.exact, prob.t0, prob.tf, Nfine);
[tc, yc] = eulerSolve(prob.f, prob.t0, prob.tf, prob.y0, Ncoarse);

figure('Name', sprintf('Euler Comparison - Case %d', caseId));
plot(tf, yexactFine, 'LineWidth', 1.8); hold on;
plot(tc, yc, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
grid on;
xlabel('t'); ylabel('y(t)');
legend('Exact solution (fine mesh)', sprintf('Euler (N = %d)', Ncoarse), 'Location', 'best');
title(sprintf('Forward Euler: Case %d', caseId));

% 3) Global error vs h
Ns = 2.^(1:4);
hs = zeros(size(Ns));
errs = zeros(size(Ns));
for k = 1:numel(Ns)
    N = Ns(k);
    [t, y] = eulerSolve(prob.f, prob.t0, prob.tf, prob.y0, N);
    ye = prob.exact(t);
    hs(k) = (prob.tf - prob.t0) / N;
    errs(k) = max(abs(ye - y));
end

figure('Name', sprintf('Euler Convergence - Case %d', caseId));
loglog(hs, errs, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6);
grid on;
xlabel('h'); ylabel('Global error ||e||_\infty');
title(sprintf('Forward Euler Convergence: Case %d', caseId));

end

function prob = getProblem(caseId)
switch caseId
    case 1
        prob.t0 = 0; prob.tf = 1; prob.y0 = 0;
        prob.f = @(t,y) t.*exp(3*t) - 2*y;
        prob.exact = @(t) (1/5).*t.*exp(3*t) - (1/25).*exp(3*t) + (1/25).*exp(-2*t);
    case 2
        prob.t0 = 1; prob.tf = 2; prob.y0 = 2;
        prob.f = @(t,y) 1 + y./t;
        prob.exact = @(t) t.*log(t) + 2*t;
    case 3
        prob.t0 = 2; prob.tf = 3; prob.y0 = 1;
        prob.f = @(t,y) 1 + (t - y).^2;
        prob.exact = @(t) t + 1./(1 - t);
    case 4
        prob.t0 = 0; prob.tf = 1; prob.y0 = 1;
        prob.f = @(t,y) cos(2*t) + sin(3*t);
        prob.exact = @(t) 0.5*sin(2*t) - (1/3)*cos(3*t) + 4/3;
    otherwise
        error('caseId must be 1, 2, 3, or 4.');
end
end

function [t, y] = eulerSolve(f, t0, tf, y0, N)
h = (tf - t0) / N;
t = linspace(t0, tf, N+1).';
y = zeros(N+1, 1);
y(1) = y0;
for n = 1:N
    y(n+1) = y(n) + h * f(t(n), y(n));
end
end

function [t, y] = exactOnGrid(exactFun, t0, tf, N)
t = linspace(t0, tf, N+1).';
y = exactFun(t);
end

function printIterations(t, y, exactFun, header)
fprintf('\n%s\n', header);
fprintf('%5s %14s %20s %20s %20s\n', 'n', 't_n', 'Euler y_n', 'Exact y(t_n)', 'Abs error');
for n = 1:numel(t)
    ye = exactFun(t(n));
    fprintf('%5d %14.8f %20.12f %20.12f %20.12e\n', n-1, t(n), y(n), ye, abs(ye - y(n)));
end
end
