function x = solvechol(A,b)

%  cholesky decomposition for symetric positive definite (SPD)

% --- Check if A_spd is symmetric ---
if ~issymmetric(A)
    error('Given Matrix A is not symmetric. Cholesky decomposition requires SPD matrix.');
end

% --- Check positive definiteness using eigenvalues ---
eig_vals = eig(A);
if any(eig_vals <= 0)
    error('Given Matrix A is not positive definite. Cholesky decomposition requires SPD matrix.');
end


L_chol = chol(A,'lower');

% Forward substitution: solve L * y = b
y = L_chol \ b;

% Back substitution: solve L' * x = y
x = L_chol' \ y;


end