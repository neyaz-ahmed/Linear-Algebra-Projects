%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB function to perform Householder replection for a matrix
function [Q, H] = householder(A)
    n = size(A, 1); 
    H = A; 
    Q = eye(n); 
    for k = 1:(n-2)
          x = H(k+1:n, k);        
        % Construct Householder vector v to zero entries of x below the first element
        norm_x = norm(x); 
        v = x + sign(x(1)) * norm_x * [1; zeros(length(x)-1, 1)]; 
        v_norm = norm(v);
        if v_norm == 0
            continue; 
        end
        u = v / v_norm;         
        w = [zeros(k, 1); u];
        Q_k = eye(n) - 2 * (w * w');
        H = Q_k' * H * Q_k;
        Q = Q * Q_k;
    end
end
