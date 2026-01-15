
clear all;
clc;

n = 10; 
% A = spdiags([-ones(n,1) 3*ones(n,1) -ones(n,1)], -1:1, n, n);
%b = randn(n,1); 
% b = b/norm(b);
%%
A = zeros(n, n);
%A= rand(n,n);
for i = 1:n
    A(i,i) = 40; 
    if i > 1
        A(i,i-1) = -1; 
        A(i-1,i) = -1; 
    end
end
[n, m] = size(A);
maxIter = 500;         
tol = 1e-1;          
%% QR Algorithm Without Shifts
Ak = A;
for k = 1:maxIter
    [Q, R] = qr(Ak);
    Ak = R * Q;
    off_diag_norm = norm(Ak - diag(diag(Ak)), 'fro');
    if off_diag_norm < tol
        fprintf('Converged without shifts after %d iterations.\n', k);
        break;
    end
end
eig1=sort(diag(full(Ak)));
%% QR Algorithm With Shifts
Ak2 = A;
for k = 1:maxIter
    mu = 6; 
    [Q, R] = qr(Ak2 - mu * speye(n));
    Ak2 = R * Q + mu * speye(n);
    off_diag_norm = norm(Ak2 - diag(diag(Ak2)), 'fro');
    if off_diag_norm < tol
        fprintf('Converged with shifts after %d iterations.\n', k);
        break;
    end
end
eig2=sort(diag(full(Ak2)));
%%  QR iteration on Hessenberg matrix
[eig3, eigvecs] = hessQR(A);
eig4=eig(A);
fprintf('Comparison with original eigen values');
disp(norm(eig1-eig4))
disp(norm(eig2-eig4))
disp(norm(eig3-eig4))



