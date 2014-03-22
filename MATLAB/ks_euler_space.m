
function ks_euler_space
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time

% Initial values
L = 32*pi;
M_sol = 2^10;
N = 10000;
k = 0.01;

U_sol = ETD_KT(M_sol, k);
error_norm = zeros(6,1);
h_p = zeros(6,1);


for j = 3:8

    M = 2^j;
    h = L/(M-1);

    % x = h * (1:M)
    % x = 32*pi/(M-1) * (1:M)
    x = (32*pi)*(1:M)/(M);

    % Boundary conditions
    % u(x,0) = f(x)
    % u(0,t) = u(L,t)

    f = @(x) cos(x/16) .* (1 + sin(x/16));

    % Creating the U-matrix and inserting boundary condition
    U = zeros(M, N);
    U(:,1) = f(x');

    % Construction of matrices for U(n+1) = (I - A - B)*U(n) - D*(U(n)).^2 
    e = ones(M,1);
    diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
    A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

    diagVecB = [-M+1 -1:1 M-1];
    B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

    diagVecD = [-M+1 -1 1 M-1];
    D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);

    % Iteration over time
    for n = 1:N
        U(:,n+1) = (eye(M)-A-B)*U(:,n) - D*(U(:,n).^2);
    end
    
    error = zeros(1,M);
    for i = 0:M-1
        error(i+1) = U(i+1,N) - U_sol((M_sol-1)/(M-1)*i+1,N);
    end
    
    error_norm(j-2) = norm(error, Inf);
    h_p(j-2) = h;
    j
end


    figure
    %contourf(t, x, U)
    loglog(h_p, error_norm, 'r');
    
end

