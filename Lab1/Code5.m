% Data
l  = [2; 4; 6];
F  = [7.0; 9.4; 12.3];
lu = 5.3;

% Design column
A = l - lu;          % 3x1 vector

% Normal system
ATA = A' * A;        % scalar (1x1)
ATF = A' * F;        % scalar (1x1)

% Cholesky: MATLAB's chol by default returns upper R such that R'*R = ATA
R = chol(ATA);       % R is 1x1 here (scalar) and R'*R = ATA

% Solve R'*y = ATF  (forward)
y = R' \ ATF;

% Solve R*k = y     (backward)
k = R \ y;

fprintf('ATA = %.6f, ATF = %.6f\n', ATA, ATF);
fprintf('Cholesky factor R = %.6f\n', R);
fprintf('k (from Cholesky) = %.6f\n', k);
