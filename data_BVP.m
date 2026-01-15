
function [A,b,x] = data_BVP()

% Parameters
n = 1000;                % number of grid points (excluding boundaries)
alpha = 0;             % boundary condition at x=0
beta  = 2;             % boundary condition at x=1
x0 = 0; x1 = 3;        % domain [x0, x1]
h = (x1 - x0) / n;     % step size

% Grid points (excluding boundaries)
x = linspace(x0 + h, x1 - h, n-1)';

% Example coefficient and source term
c = 5 * ones(n-1, 1);          % coefficient c_i
f = sin(pi * x);               % right-hand side f_i

% Construct tridiagonal matrix A
main_diag = (2 + h^2 * c) / h^2;   % main diagonal
off_diag  = -ones(n-2, 1) / h^2;   % sub- and super-diagonal

A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

% Construct vector b
b = f;
b(1)   = b(1)   + alpha / h^2;    % adjust for left boundary
b(end) = b(end) + beta  / h^2;    % adjust for right boundary

end