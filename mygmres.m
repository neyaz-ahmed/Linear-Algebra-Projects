function [x, res] = mygmres(A, b, x0, tol)
    m = length(b);
    x = x0;
    r = b - A*x; 
    b_norm = norm(b);
    r_norm = norm(r);
    res = [r_norm / b_norm]; 
    max_iter =500;
    Q = zeros(m, max_iter+1);
    Q(:,1) = r / r_norm;   % First basis vector
    H = zeros(max_iter+1, max_iter); % Hessenberg matrix
    beta = zeros(max_iter+1, 1); % Residual vector for least squares
    beta(1) = r_norm;
    for k = 1:max_iter
        % Arnoldi iteration
        [H(1:k+1, k), Q(:, k+1)] = arnoldi(A, Q, k);        
        % Solve least squares problem: minimize ||beta*e1 - H*y||
        e1 = [1; zeros(k, 1)];
        y = H(1:k+1, 1:k) \ (beta(1) * e1);        
        % Compute approximate solution
        x = x0 + Q(:, 1:k) * y;        
        % Compute residual norm
        r = b - A*x;
        res_norm = norm(r);
        res = [res; res_norm / b_norm];        
        % Check convergence
        if res_norm < tol
            fprintf('GMRES converged after %d iterations with residual norm %e\n', k, res_norm);
            break;
        end        
        % Check for breakdown
        if H(k+1, k) < 1e-12
            fprintf('Breakdown in Arnoldi iteration at step %d\n', k);
            break;
        end
    end    
end
