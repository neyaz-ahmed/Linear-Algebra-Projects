function [lambda, X, iter] = inversepower(A)
%=========================================================
% Inverse Iteration with Rayleigh Quotient
%---------------------------------------------------------
%
% INPUT:  A  : (n x n) matrix (can be unsymmetric)
%
% OUTPUT:
%   lambda   : vector of computed eigenvalues
%   X        : matrix of corresponding eigenvectors
%   iter     : number of iterations for each eigenpair
%%---------------------------------------------------------
n = size(A, 1);
lambda = zeros(n, 1);
X = zeros(n, n);
iter = zeros(n, 1);
max_iter = 1000;
tol = 1e-10;
for i = 1:n
    x = randn(n, 1);
    if i > 1
        x = x - X(:, 1:i-1) * (X(:, 1:i-1)' * x);
    end
    x = x / norm(x);
    mu = (x' * A * x) / (x' * x);
    for k = 1:500
        y = (A - mu * eye(n)) \ x;
        x_new = y / norm(y);
        lambda_i = (x_new' * A * x_new) / (x_new' * x_new);
        r = norm(A * x_new - lambda_i * x_new);
        if r < tol
            lambda(i) = lambda_i;
            X(:, i) = x_new;
            iter(i) = k;
            break;
        end
        x = x_new;
        mu = lambda_i;
    end
    if iter(i) == 0
        warning('Eigenpair %d did not converge.', i);
    end
end
end
