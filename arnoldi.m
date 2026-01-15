%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB function to perform Lanczos algorithm for a matrix
function [V, H] = arnoldi(A, b, m)
n = size(A, 1);
epsilon = 1e-12; 
V = zeros(n, m+1); 
H = zeros(m+1, m); 
k = m; 
V(:, 1) = b / norm(b);     
    for i = 1:m
        w = A * V(:, i);               
   % Reorthogonalize w against all previous vectors (Gram-Schmidt)
        for j = 1:i
            h = V(:, j)' * w;
            w = w - h * V(:, j);
            H(j, i) = H(j, i) + h; 
        end        
        H(i+1, i) = norm(w);              
        if H(i+1, i) < epsilon
            fprintf('Arnoldi stopped at iteration %d: h_%d,%d < epsilon.\n', i, i+1, i);
            k = i; 
            break;
        end
        V(:, i+1) = w / H(i+1, i);
    end
    V = V(:, 1:k+1); 
    H = H(1:k+1, 1:k);   
end

