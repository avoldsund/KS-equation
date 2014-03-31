function [error_norm] = IMEX_conv_time()
% Function to check convergence in time for the implicit-explicit scheme

global Ms hs;
% Reference solution

M = 2^8;
length = 32*pi;
h = length/M;
hs = h;
Ms = M;
x = h * (0:M-1);

time = 10;
N_sol = 2^14;
k_sol = time/N_sol;
U_sol = ref_sol(k_sol, time, x);

low = 7;
high = 12;

error_norm = zeros(high-low+1, 1);
kVec = zeros(high-low+1, 1);

for j = low:high
    N = 2^j;
    k = time/N;
    U = zeros(M, N);
    f = @(x) cos(x/16).*(1+sin(x/16));
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
    error = U(:,end) - U_sol(:,end);
    error_norm(j-low+1, 1) = norm(error,Inf);
    kVec(j-low+1, 1) = k; 
end

% Loglog-plot of the error
figure
loglog(kVec, error_norm, 'r', kVec, kVec, 'g')
xlabel('Log of the time steps k')
ylabel('Log of the error norm')
title('Loglog-plot of the error norm')
legend('Infinity-norm','Slope 1','location','SouthEast')

end