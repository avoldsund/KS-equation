tic
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time

% Initial values
length = 32*pi;
M = 128;
h = length/(M-1);

time = 30;
N = 30000;
k = time/N;

% x = h * (1:M)
% x = 32*pi/(M-1) * (1:M)
x = (32*pi)*(1:M)/(M);
t = k * (0:N);

% Boundary conditions
% u(x,0) = f(x)
% u(0,t) = u(L,t)
 
f = @(x) cos(x/16) .* (1 + sin(x/16));

% Creating the U-matrix and inserting boundary condition
U = zeros(M, N);
U(:,1) = f(x');

% Construction of the A, B and D matrix
e = ones(M,1);
diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

diagVecB = [-M+1 -1:1 M-1];
B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

diagVecD = [-M+1 -1 1 M-1];
D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);

% Running 
for i = 1:N
    U(:,i+1) = U(:,i) - A*U(:,i) - B*U(:,i) - D*(U(:,i).^2);
end

figure
%contourf(t, x, U)
contourf(U')
toc