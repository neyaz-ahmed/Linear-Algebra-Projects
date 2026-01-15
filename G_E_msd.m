clear; clc; close all;

%% System matrices (second-order form)
m1 = 1;  m2 = 1;
k1 = 3;  k2 = 2;  k3 = 3;
c1 = 0.2; c2 = 0.2;

M = [m1 0;
      0 m2];                           
C = [c1+c2  -c2;  -c2  c2];                    
K = [k1+k2  -k2;  -k2  k2+k3];                  

%% Build the state-space matrix A
n = 2;
A = [ zeros(n)    eye(n);
     -M\K       -M\C ];

%% Compute eigenvalues and eigenvectors
[V, Lambda] = eig(A,'vector');                 
[eigval, idx] = sort(Lambda,'ComparisonMethod','real');
V = V(:,idx);                    
%% Initial conditions 
x0  = [0.1; -0.05];     
v0  = [0; 0];           
z0  = [x0; v0];         

c = V \ z0;   

tspan = 0:0.1:5;
%% Exact analytical solution
Z_exact = zeros(length(tspan),4);
for i = 1:length(tspan)
    t = tspan(i);
    Z_exact(i,:) = real(V * (c .* exp(eigval*t)) );  
end
x1_exact = Z_exact(:,1);
x2_exact = Z_exact(:,3);
%% Numerical integration with ode45 (for comparison)
ode_fun = @(t,z) A*z;
[~, Z_num] = ode45(ode_fun, tspan, z0);
x1_num = Z_num(:,1);
x2_num = Z_num(:,3);
%% Plot results
figure(5);
plot(tspan, x1_exact, 'b-', 'LineWidth', 1.5, 'DisplayName','Exact (eig)');
hold on;
plot(tspan, x1_num, 'r--', 'LineWidth', 1.5, 'DisplayName','ode45');
ylabel('Displacement'); 
grid on;

figure(6);
plot(tspan, x2_exact, 'b-', 'LineWidth', 2);
hold on;
plot(tspan, x2_num, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Velocity'); 
grid on;
legend('Exact (eigen)','ode45','Location','northeast');

