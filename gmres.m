
function x_m = gmres(A,b,x0)

%% Define system Ax = b
% n = 50;                                  % dimension of A
% A = randn(n);                           % random (non-symmetric) system
% x_exact = randn(n,1);                         % right-hand side
% b=A*x_exact;
% x0 = zeros(n,1);                        % initial guess
% tol = 1e-6;                             % tolerance
max_iter = size(A,1);                           % maximum dimension of Krylov subspace (restart value m)

%% Step 1: Compute initial residual
r0 = b - A*x0;
beta = norm(r0,2);
v1 = r0 / beta;
V = v1;                                 % first basis vector
H = zeros(max_iter+1, max_iter);        % Hessenberg matrix

%% Arnoldi Process
for j = 1:max_iter
    w = A * V(:,j);
    % Orthogonalize against existing basis vectors
    for i = 1:j
        H(i,j) = V(:,i)' * w;
        w = w - H(i,j) * V(:,i);
    end
    
    H(j+1,j) = norm(w,2);
    
    if H(j+1,j) == 0
        fprintf('Krylov subspace became invariant at step %d.\n', j);
        break;
    end    
    V = [V, w / H(j+1,j)]; 
end

%% Solve least-squares problem: min ||beta*e1 - H*y||
e1 = zeros(j+1,1);
e1(1) = 1;
rhs = beta * e1;
y = H(1:j+1,1:j) \ rhs;  

%% Compute approximate solution
x_m = x0 + V(:,1:j) * y;

%% Compute residual and check convergence
% residual = norm(x_exact - x_m,2) / norm(b,2);
% fprintf('Final relative residual: %.2e\n', residual);
% 
% if residual < tol
%     fprintf('GMRES converged in %d iterations.\n', j);
% else
%     fprintf('GMRES did not converge within %d iterations.\n', j);
% end

end