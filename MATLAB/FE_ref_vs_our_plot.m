function [U k] = FE_ref_vs_our_plot(length, Ms, time, N)
% ks_euler_anders(length, M, time, N)
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
global Ms hs
% Initial values
hs = length/Ms;
k = time/N;

x = hs * (0:Ms-1);
t = k * (0:N-1);

% Boundary conditions
% u(x,0) = f(x)
% u(0,t) = u(L,t)

f = @(x) cos(x/16) .* (1 + sin(x/16));

% Creating the U-Matrix and inserting boundary condition
U = zeros(Ms, N);
U(:,1) = f(x');

% Construction of the Axxxx-, Axx- and B-matrix
e = ones(Ms,1);
diagVecAxxxx = [-Ms+1 -Ms+2 -2:2 Ms-2 Ms-1];
Axxxx = (k/(hs^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecAxxxx, Ms, Ms);

diagVecAxx = [-Ms+1 -1:1 Ms-1];
Axx = (k/(hs^2)) * spdiags([e e -2*e e e], diagVecAxx, Ms, Ms);

diagVecB = [-Ms+1 -1 1 Ms-1];
B = (k/(4*hs)) * spdiags([1*e -1*e 1*e -1*e], diagVecB, Ms, Ms);

% For-loop to calculate the next time step
for i = 1:N
    U(:,i+1) = U(:,i) - Axxxx*U(:,i) - Axx*U(:,i) - B*(U(:,i).^2);
end

U_sol = ref_sol_global(k,time,x);


figure
hold on

plot(x,U_sol(:,end),'r')
plot(x,U(:,end))
end