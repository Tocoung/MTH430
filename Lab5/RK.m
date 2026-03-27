function rk_method_ivp(caseId, scheme)
%RK_METHOD_IVP  Runge-Kutta solvers for the 4 given IVPs.
%   rk_method_ivp(caseId, scheme)
%   scheme = 'midpoint', 'modifiedeuler', or 'rk4'.
%   Produces iteration printout, comparison plot, and convergence plot.

if nargin < 1, caseId = 1; end
if nargin < 2, scheme = 'rk4'; end
scheme = lower(char(scheme));
allowed = {'midpoint', 'modifiedeuler', 'rk4'};
ok = any(strcmp(scheme, allowed));
if ~ok
    error('scheme must be midpoint, modifiedeuler, or rk4.');
end

prob = getProblem(caseId);

Ncoarse = 2^3;
[tc, yc] = rkSolve(prob.f, prob.t0, prob.tf, prob.y0, Ncoarse, scheme);
printIterations(tc, yc, prob.exact, sprintf('RK (%s) | Case %d | N = %d', scheme, caseId, Ncoarse));

Nfine = 2^9;
[tf, yexactFine] = exactOnGrid(prob.exact, prob.t0, prob.tf, Nfine);
[tc, yc] = rkSolve(prob.f, prob.t0, prob.tf, prob.y0, Ncoarse, scheme);

figure('Name', sprintf('RK-%s Comparison - Case %d', scheme, caseId));
plot(tf, yexactFine, 'LineWidth', 1.8); hold on;
plot(tc, yc, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
grid on;
xlabel('t'); ylabel('y(t)');
legend('Exact solution (fine mesh)', sprintf('RK %s (N = %d)', scheme, Ncoarse), 'Location', 'best');
title(sprintf('Runge-Kutta %s: Case %d', scheme, caseId));

Ns = 2.^(1:4);
hs = zeros(size(Ns));
errs = zeros(size(Ns));
for k = 1:numel(Ns)
    N = Ns(k);
    [t, y] = rkSolve(prob.f, prob.t0, prob.tf, prob.y0, N, scheme);
    ye = prob.exact(t);
    hs(k) = (prob.tf - prob.t0) / N;
    errs(k) = max(abs(ye - y));
end

figure('Name', sprintf('RK-%s Convergence - Case %d', scheme, caseId));
loglog(hs, errs, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6);
grid on;
xlabel('h'); ylabel('Global error ||e||_\infty');
title(sprintf('Runge-Kutta %s Convergence: Case %d', scheme, caseId));

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

function [t, y] = rkSolve(f, t0, tf, y0, N, scheme)
h = (tf - t0) / N;
t = linspace(t0, tf, N+1).';
y = zeros(N+1, 1);
y(1) = y0;
for n = 1:N
    tn = t(n);
    yn = y(n);
    switch scheme
        case 'midpoint'
            k1 = f(tn, yn);
            k2 = f(tn + h/2, yn + (h/2)*k1);
            y(n+1) = yn + h*k2;
        case 'modifiedeuler'
            k1 = f(tn, yn);
            k2 = f(tn + h/2, yn + (h/2)*k1);
            y(n+1) = yn + (h/2)*(k1 + k2);
        case 'rk4'
            k1 = f(tn, yn);
            k2 = f(tn + h/2, yn + (h/2)*k1);
            k3 = f(tn + h/2, yn + (h/2)*k2);
            k4 = f(tn + h, yn + h*k3);
            y(n+1) = yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        otherwise
            error('Unknown scheme.');
    end
end
end

function [t, y] = exactOnGrid(exactFun, t0, tf, N)
t = linspace(t0, tf, N+1).';
y = exactFun(t);
end

function printIterations(t, y, exactFun, header)
fprintf('\n%s\n', header);
fprintf('%5s %14s %20s %20s %20s\n', 'n', 't_n', 'RK y_n', 'Exact y(t_n)', 'Abs error');
for n = 1:numel(t)
    ye = exactFun(t(n));
    fprintf('%5d %14.8f %20.12f %20.12f %20.12e\n', n-1, t(n), y(n), ye, abs(ye - y(n)));
end
end
