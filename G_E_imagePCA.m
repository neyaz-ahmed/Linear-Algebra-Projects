%%    image_compression

clear; clc; 

image_path = 'Frida.jpg'; 
A = imread(image_path);   
A = double(A) / 255; 
[m, n] = size(A);

[U, S, V] = svd(A, 'econ'); 

reduced = 100;
k = min(reduced, min(m, n)); 
A_approx = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
% A_approx = U * S* V';

compression_ratio = k * (m + n + 1) / (m * n);

disp('Compression Ratio:');
fprintf('%.4f (lower is better)\n', compression_ratio);

% Display results
figure(1);
subplot(1, 2, 1);
imagesc(A);
colormap gray;
axis equal; axis tight;
title('Original Image', 'FontSize', 12);
colorbar;
% Approximated image
subplot(1, 2, 2);
imagesc(A_approx);
colormap gray;
axis equal; axis tight;
title(['Rank-', num2str(k), ' Approximation'], 'FontSize', 12);
colorbar;

reconstruction_error = norm(A - A_approx, 'fro') / norm(A, 'fro');
disp('Relative Reconstruction Error:');
fprintf('%.4e\n', reconstruction_error);
