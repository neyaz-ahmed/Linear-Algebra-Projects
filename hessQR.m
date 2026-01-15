function [eigvals, eigvecs] = hessQR(A)
% eigenQR computes the eigenvalues of a matrix using the QR algorithm
%
   % Step 1: Reduce to Hessenberg form (speeds up QR iterations)
     [Qh, H] = hess(full(A)); 
     %[Qh, H] = myHessenberg(A);   % hess(A)
    Q_total = Qh; % Initialize eigenvector matrix
    % Parameters
    maxIter = 1000;      % Maximum iterations
    tol = 1e-12;         % Convergence tolerance
    % Step 2: QR iteration loop
    for k = 1:maxIter
        [Q, R] = qr(H);  % QR decomposition
        H = R * Q;       % Update Hessenberg matrix
        Q_total = Q_total * Q; % Accumulate orthogonal matrices for eigenvectors
        % Step 3: Check convergence (subdiagonal entries close to zero)
        if norm(tril(H, -1), 'fro') < tol
            fprintf('QR algorithm converged after %d iterations.\n', k);
            break;
        end
    end
    
    % Step 4: Extract eigenvalues from diagonal
    eigvals = diag(H);
    eigvecs = Q_total; % Eigenvectors are columns of Q_total
    % If not converged, warn the user
    if k == maxIter
        warning('QR algorithm did not fully converge within %d iterations.', maxIter);
    end
end
