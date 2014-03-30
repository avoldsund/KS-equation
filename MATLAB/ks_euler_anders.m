function [U h x k] = ks_euler_anders(length, M, time, N)
% ks_euler_anders(length, M, time, N)
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time

% Initial values
h = length/M;
k = time/N;

x = h * (0:M-1);
t = k * (0:N-1);

% Boundary conditions
% u(x,0) = f(x)
% u(0,t) = u(L,t)

f = @(x) cos(x/16) .* (1 + sin(x/16));

% Creating the U-matrix and inserting boundary condition
U = zeros(M, N);
U(:,1) = f(x');

% Construction of the Axxxx-, Axx- and B-matrix
e = ones(M,1);
diagVecAxxxx = [-M+1 -M+2 -2:2 M-2 M-1];
Axxxx = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecAxxxx, M, M);

diagVecAxx = [-M+1 -1:1 M-1];
Axx = (k/(h^2)) * spdiags([e e -2*e e e], diagVecAxx, M, M);

diagVecB = [-M+1 -1 1 M-1];
B = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecB, M, M);

% For-loop to calculate the next time step
for i = 1:N-1
    U(:,i+1) = U(:,i) - Axxxx*U(:,i) - Axx*U(:,i) - B*(U(:,i).^2);
end
% 
% [V, xa, ya] = compress(t, x, U, 1000, 1000);
% 
% figure
% mesh(xa,ya,V')
% %view(230,45)
% view(340,60)
% %view(90,270) Contour plot of the solution
% title('Numerical solution of the KS-eq: M = 128, N = 2^{16}')
% ylabel('time')
% xlabel('space [0, 32*pi]')

% figure
% surf(t,x,U)
% title('Numerical solution for KS-equation')
% xlabel('time')
% ylabel('space')
% view(90,270)
% shading FLAT
% %contourf(t,x,U)
end
