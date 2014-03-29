global Ms hs

Ms = 2^10;
length = 32*pi;
hs = length/Ms;
x = hs * (0:Ms-1);
r = k/h^2;

time = 10;
N_sol = 2^20;
k_sol = time/N_sol;
U_sol = ref_sol(k_sol, time, x);