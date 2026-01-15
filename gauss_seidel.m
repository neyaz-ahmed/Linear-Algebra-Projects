 function x = gauss_seidel(A, b)
    n = length(b);
    x = rand(n, 1);           % Initial guess;
    max_iter= 500;
    for k = 1:max_iter
        x_new = x;
        for i = 1:n
            sum_ax = 0;
            for j = 1:n
                if j < i
                    sum_ax = sum_ax + A(i,j) * x_new(j);
                elseif j > i
                    sum_ax = sum_ax + A(i,j) * x(j);
                end
            end
            %%  EDIT YOU
            x_new(i) = (b(i) - sum_ax) / A(i,i);
        end
        if norm(x_new - x, inf) < 1e-6
            x = x_new;
            fprintf('Gauss-Seidel converged in %d iterations\n', k);
            return;
        end
        x = x_new;
    end
    warning('Gauss-Seidel method did not converge within %d iterations', max_iter);
end