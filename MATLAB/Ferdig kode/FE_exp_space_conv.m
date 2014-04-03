function [error_norm h_p] = FE_exp_space_conv

% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
% Reference solution with ode15

    % Boundary conditions
    % u(x,0) = f(x)
    % u(0,t) = u(L,t)

f = @(x) cos(x/16).*(1+sin(x/16));

global Ms hs
L = 32*pi;
Ms = 2^11;
hs = L/Ms;
k = 0.0001;
T = 10;
N = T/k;
xs = 0:hs:L-hs;

% Reference solution

yy = ref_sol(k,T,xs);

disp('done')

% Parameters for for loop
min = 6;
max = 9;

num = max-min+1;
error_norm = zeros(num,1);
h_p = zeros(num,1);


for j = min:max
    j
    M = 2^j;
    h = L/M;
    x = 0:h:L-h;
    
    % Creating the U-matrix and inserting boundary condition
    U = zeros(M, N);
    size(U)
    size(x)
    U(:,1) = f(x');

    % Construction of matrices for U(n+1) = (I - A - B)*U(n) - D*(U(n)).^2 
    A = k/(h^2)*second_order_matrix(M);
    B = k/(h^4)*second_order_matrix(M)*second_order_matrix(M);
    D = k/(4*h)*first_order_central_matrix(M);
    
    % Iteration over time
    for n = 1:N
        U(:,n+1) = (eye(M)-A-B)*U(:,n) - D*(U(:,n).^2);
    end
    
    error = (yy(1:Ms/M:Ms,N)-U(:,N));
    error_norm(j-min+1) = norm(error, Inf);
    h_p(j-min+1) = h;
end

figure
    loglog(h_p, error_norm, 'ro-', h_p, h_p, 'b', h_p, h_p.^2, 'g');
    
end
