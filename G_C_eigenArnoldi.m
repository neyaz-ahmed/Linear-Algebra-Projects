% % ------------------------------------------------------
% Test Arnoldi and Lanczos methods for large-scale eigenvalue computation
% ------------------------------------------------------
clc; clear; close all;

%%
 n = 10; m = 10; epsilon = 1e-12;
 A = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], -1:1, n, n);
b = randn(n,1); b = b/norm(b);
%% Arnoldi iteration (non-symmetric case)
disp('=== Arnoldi Iteration ===');
[V_arn, H_arn] = arnoldi(A, b, m);

[vecs_H, D] = eig(H_arn(1:m,1:m));
disp('Ritz eigenvalues (Arnoldi approximation):');
ritz_vals_arn = sort(diag(D), 'descend'); % Ritz values
disp(ritz_vals_arn);  % Display a few
% Approximate eigenvectors from Arnoldi
approx_vecs_arn = V_arn(:,1:m) * vecs_H;

%% Lanczos iteration (symmetric case)
disp('=== Lanczos Iteration ===');
[V_lan, T_lan] = lanczos(A, b, m);

 [vecs_T, D] = eig(T_lan(1:m,1:m));
disp('Ritz eigenvalues (Lanczos approximation):');
ritz_vals_lan = sort(diag(D), 'descend'); % Ritz values
disp(ritz_vals_lan);  % Display a few
% Approximate eigenvectors from Lanczos
approx_vecs_lan = V_lan(:,1:m) * vecs_T;

%% Compare with MATLAB eigs (reference for large sparse matrices)
disp('=== MATLAB eigs (reference) ===');
ref_vals = (eigs(A)); %, 6, 'largestabs');
disp('Dominant eigenvalues from eigs:');
disp(ref_vals);

%% Check orthogonality of Arnoldi/Lanczos bases
disp('Check orthogonality of Arnoldi basis (Q''*Q):');
disp(norm(V_arn(:,1:m)' * V_arn(:,1:m)));

disp('Check orthogonality of Lanczos basis (Q''*Q):');
disp(norm(V_lan(:,1:m)' * V_lan(:,1:m)));
