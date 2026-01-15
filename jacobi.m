function x = jacobi(A, b)
    n = length(b);
    x = rand(n, 1);        % Initial guess x0;
    max_iter=1000;
    for k = 1:max_iter
        x_new = zeros(n, 1);
        for i = 1:n
            sum_ax = 0;
            for j = 1:n
                if j ~= i
                    sum_ax = sum_ax + A(i,j) * x(j);
                end
            end
            x_new(i) = (b(i) - sum_ax) / A(i,i);
        end
        if norm(x_new - x) < 1e-6
            x = x_new;
            fprintf('Jacobi converged in %d iterations\n', k);
            return;
        end
        x = x_new;
    end
warning('Jacobi method did not converge within %d iterations', max_iter);
end