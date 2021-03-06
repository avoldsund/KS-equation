function [error_norm error_norm_unscaled] = FE_exp_time_conv(U_sol)
% Prove convergence in time for the KS-equation
global Ms hs;
% Finding the reference solution
M = 2^9;
length = 32*pi;
h = length/M;
x = h * (0:M-1);
hs = h;
Ms = M;

time = 10;
N_sol = 2^14;
k_sol = time/N_sol;
U_sol = ref_sol(k_sol, time, x);

low = 7;
high = 12;

error_norm = zeros(high-low+1, 1);
error_norm_unscaled = zeros(high-low+1, 1);
kVec = zeros(high-low+1, 1);

for j = low:high
    N = 2^j;
    k = time/N;
%     [U k] = ks_euler_anders(length, Ms, time, N);
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
        U(:,n+1) = (eye(M)-A-B)*U(:,n) - D*(U(:,n).^2);
%         U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    end
    
    
    kVec(j-low+1, 1) = k;   
    
    error = U(:,end) - U_sol(:,end);
%     error_norm_unscaled(j-low+1, 1) = norm(error,Inf);
    error_norm(j-low+1, 1) = norm(error,Inf);
    
end

loglog(kVec, error_norm, 'r', kVec, kVec, 'g')

end