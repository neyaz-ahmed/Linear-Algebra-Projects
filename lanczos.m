%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB function to perform Lanczos algorithm for a matrix
function [V, T] = lanczos(A, b, m)
    n = size(A, 1);
    epsilon = 1e-12; 
    V = zeros(n, m); 
    alpha = zeros(m, 1); 
    beta = zeros(m, 1); 
    v0 = zeros(n, 1); 
    v1 = b / norm(b); 
    V(:, 1) = v1;    
    k = m; 
    for j = 1:m
        w = A * v1- beta(j) * v0;        
        alpha(j) = v1' * w;        
        w = w - alpha(j) * v1;        
   % Reorthogonalize w against all previous vectors (Gram-Schmidt)
        for i = 1:j
            w = w - (V(:, i)' * w) * V(:, i);
        end
        beta(j+1) = norm(w);
        if beta(j+1) < epsilon
            fprintf('Lanczos stopped at iteration %d: T_%d,%d < epsilon.\n', j, j+1, j);
            k = j; 
            break;
        end
        v0 = v1;
        v1 = w / beta(j+1);
        if j < m
            V(:, j+1) = v1;
        end
    end
    T = diag(alpha(1:k));
    if k > 1
        T = T + diag(beta(2:k), 1) + diag(beta(2:k), -1);
    end          
end

