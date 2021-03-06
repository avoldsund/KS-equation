function [error_norm h_p] = IMEX_conv_space()
% Function to check convergence in space for the implicit-explicit scheme

f = @(x) cos(x/16).*(1+sin(x/16));

global Ms hs
L = 32*pi;
Ms = 2^13;
hs = L/Ms;

% Reference solution
k = 0.0005;
xx = 0:hs:L-hs;
N = 100000;
T = N*k;
yy = ref_sol(k,T,xx);

% M values increasing by a factor 2: 2^j
min = 7;
max = 10;
num = max-min+1;
error_norm = zeros(num,1);
h_p = zeros(num,1);

for j = min:max

    M = 2^j;
    h = L/M;
    x = 0:h:L-h;

    % Creating the U-matrix and inserting boundary condition
    U = zeros(M, N);
    U(:,1) = f(x');

    % Construction of matrices for U(n+1) = (I - A - B)*U(n) - D*(U(n)).^2 

    A = k/(h^2)*second_order_matrix(M);
    B = k/(h^4)*second_order_matrix(M)*second_order_matrix(M);
    D = k/(4*h)*first_order_central_matrix(M);
    
    F = (speye(M)+A/2+B/2);
    G = (speye(M)-A/2-B/2);
    
    % Iteration over time
    for n = 1:N
        U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    end

    % Finding the norm of the error
    error = (yy(1:Ms/M:Ms,N)-U(:,N+1));
    error_norm(j-min+1) = norm(error, Inf);
    h_p(j-min+1) = h;
end

% Loglog-plot of the error
figure
loglog(h_p, error_norm, 'ro-', h_p, h_p.^2, 'g');
xlabel('Space step h')
ylabel('Error norm')
title('Implicit scheme - Convergence in space')
legend('error(h)','slope 1','location','SouthEast')
    
end
