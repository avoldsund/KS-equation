function [error_norm error_norm_unscaled] = FE_exp_time_conv(U_sol)
% Prove convergence in time for the KS-equation
global Ms hs;
% Finding the reference solution
Ms = 2^9;
length = 32*pi;
hs = length/Ms;
x = hs * (0:Ms-1);

time = 10;
N_sol = 2^20;
k_sol = time/N_sol;
U_sol = ref_sol(k_sol, time, x);

low = 16;
high = 20;

error_norm = zeros(high-low+1, 1);
error_norm_unscaled = zeros(high-low+1, 1);
kVec = zeros(high-low+1, 1);

for j = low:high
    N = 2^j;
    [U k] = ks_euler_anders(length, Ms, time, N);
    kVec(j-low+1, 1) = k;   
    
    error = U(:,end) - U_sol(:,end);
    error_norm_unscaled(j-low+1, 1) = norm(error,2);
    error_norm(j-low+1, 1) = sqrt(k) * norm(error,2);
    
end

loglog(kVec, error_norm, 'r', kVec, kVec, 'g')

end