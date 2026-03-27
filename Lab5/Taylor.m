function taylor_method_ivp(caseId, order)
%TAYLOR_METHOD_IVP  Taylor solver of order 2 or 4 for the 4 given IVPs.
%   taylor_method_ivp(caseId, order)
%   order = 2 or 4.
%   Produces iteration printout, comparison plot, and convergence plot.

if nargin < 1, caseId = 1; end
if nargin < 2, order = 2; end
if ~ismember(order, [2, 4])
    error('order must be 2 or 4.');
end

prob = getProblem(caseId);

Ncoarse = 2^3;
[tc, yc] = taylorSolve(caseId, prob.t0, prob.tf, prob.y0, Ncoarse, order);
printIterations(tc, yc, prob.exact, sprintf('Taylor order %d | Case %d | N = %d', order, caseId, Ncoarse));

Nfine = 2^9;
[tf, yexactFine] = exactOnGrid(prob.exact, prob.t0, prob.tf, Nfine);
[tc, yc] = taylorSolve(caseId, prob.t0, prob.tf, prob.y0, Ncoarse, order);

figure('Name', sprintf('Taylor-%d Comparison - Case %d', order, caseId));
plot(tf, yexactFine, 'LineWidth', 1.8); hold on;
plot(tc, yc, 'o-', 'LineWidth', 1.4, 'MarkerSize', 5);
grid on;
xlabel('t'); ylabel('y(t)');
legend('Exact solution (fine mesh)', sprintf('Taylor order %d (N = %d)', order, Ncoarse), 'Location', 'best');
title(sprintf('Taylor order %d: Case %d', order, caseId));

Ns = 2.^(1:4);
hs = zeros(size(Ns));
errs = zeros(size(Ns));
for k = 1:numel(Ns)
    N = Ns(k);
    [t, y] = taylorSolve(caseId, prob.t0, prob.tf, prob.y0, N, order);
    ye = prob.exact(t);
    hs(k) = (prob.tf - prob.t0) / N;
    errs(k) = max(abs(ye - y));
end

figure('Name', sprintf('Taylor-%d Convergence - Case %d', order, caseId));
loglog(hs, errs, 'o-', 'LineWidth', 1.6, 'MarkerSize', 6);
grid on;
xlabel('h'); ylabel('Global error ||e||_\infty');
title(sprintf('Taylor order %d Convergence: Case %d', order, caseId));

end

function prob = getProblem(caseId)
switch caseId
    case 1
        prob.t0 = 0; prob.tf = 1; prob.y0 = 0;
        prob.exact = @(t) (1/5).*t.*exp(3*t) - (1/25).*exp(3*t) + (1/25).*exp(-2*t);
    case 2
        prob.t0 = 1; prob.tf = 2; prob.y0 = 2;
        prob.exact = @(t) t.*log(t) + 2*t;
    case 3
        prob.t0 = 2; prob.tf = 3; prob.y0 = 1;
        prob.exact = @(t) t + 1./(1 - t);
    case 4
        prob.t0 = 0; prob.tf = 1; prob.y0 = 1;
        prob.exact = @(t) 0.5*sin(2*t) - (1/3)*cos(3*t) + 4/3;
    otherwise
        error('caseId must be 1, 2, 3, or 4.');
end
end

function [t, y] = taylorSolve(caseId, t0, tf, y0, N, order)
h = (tf - t0) / N;
t = linspace(t0, tf, N+1).';
y = zeros(N+1, 1);
y(1) = y0;
for n = 1:N
    tn = t(n);
    yn = y(n);
    [y1, y2, y3, y4] = derivatives(caseId, tn, yn);
    if order == 2
        y(n+1) = yn + h*y1 + (h^2/2)*y2;
    else
        y(n+1) = yn + h*y1 + (h^2/2)*y2 + (h^3/6)*y3 + (h^4/24)*y4;
    end
end
end

function [y1, y2, y3, y4] = derivatives(caseId, t, y)
switch caseId
    case 1
        y1 = t.*exp(3*t) - 2*y;
        y2 = (1 + 3*t).*exp(3*t) - 2*y1;
        y3 = (6 + 9*t).*exp(3*t) - 2*y2;
        y4 = (27 + 27*t).*exp(3*t) - 2*y3;
    case 2
        y1 = 1 + y./t;
        y2 = 1./t;
        y3 = -1./(t.^2);
        y4 = 2./(t.^3);
    case 3
        u = t - y;
        y1 = 1 + u.^2;
        y2 = -2*u.^3;
        y3 = 6*u.^4;
        y4 = -24*u.^5;
    case 4
        y1 = cos(2*t) + sin(3*t);
        y2 = -2*sin(2*t) + 3*cos(3*t);
        y3 = -4*cos(2*t) - 9*sin(3*t);
        y4 = 8*sin(2*t) - 27*cos(3*t);
    otherwise
        error('caseId must be 1, 2, 3, or 4.');
end
end

function [t, y] = exactOnGrid(exactFun, t0, tf, N)
t = linspace(t0, tf, N+1).';
y = exactFun(t);
end

function printIterations(t, y, exactFun, header)
fprintf('\n%s\n', header);
fprintf('%5s %14s %20s %20s %20s\n', 'n', 't_n', 'Taylor y_n', 'Exact y(t_n)', 'Abs error');
for n = 1:numel(t)
    ye = exactFun(t(n));
    fprintf('%5d %14.8f %20.12f %20.12f %20.12e\n', n-1, t(n), y(n), ye, abs(ye - y(n)));
end
end
