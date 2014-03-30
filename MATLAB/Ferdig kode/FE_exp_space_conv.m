function error_norm = FE_exp_space_conv
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
% Reference solution with ode15

    % Boundary conditions
    % u(x,0) = f(x)
    % u(0,t) = u(L,t)

f = @(x) cos(x/16).*(1+sin(x/16));

global Ms hs

% Constants
L = 32*pi;
T = 10;
% N = 2^14;
% k = T/N;
k = 10^(-4);
N = T/k;
y = 0:hs:L-hs;

% Reference solution
Ms = 2^10;
hs = L/Ms;
yy = ref_sol(k,T,y);

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
    h = L/(M);
    %k/(h^2) <= 0.0026 for convergence
    x = 0:h:L-h;

    % Creating the U-matrix and inserting boundary condition
    U = zeros(M, N);
    U(:,1) = f(x');

    % Construction of matrices for U(n+1) = (I - A - B)*U(n) - D*(U(n)).^2 

    A = k/(2*h^2)*second_order_matrix(M);
    B = k/(2*h^4)*second_order_matrix(M)*second_order_matrix(M);
    D = k/(4*h)*first_order_central_matrix(M);
    
%     F = (speye(M)+A+B);
%     G = (speye(M)-A-B);
    
    % Iteration over time
    for n = 1:N
        U(:,n+1) = (eye(M)-A-B)*U(:,n) - D*(U(:,n).^2);
%         U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    end


    % error compared to reference solution
%     error = zeros(1,M);
%     for i = 0:M-1
%         error(i+1) = U(i+1,N) - yy((Ms)/(M)*i+1,N);
%     end
    error = (yy(1:Ms/M:Ms,N)-U(:,N));
    error_norm(j-min+1) = norm(error, Inf);
    h_p(j-min+1) = h;
end

figure
    loglog(h_p, error_norm, 'ro-', h_p, h_p, 'b', h_p, h_p.^2, 'g');
    
end