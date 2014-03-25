tic
% Implementation of the Kuramoto-Sivashinsky equation
% u_t + u_xxxx + u_xx + uu_x = 0
% Central differences in space, forwards in time
% U(n+1) = (I - A - B)U(n) - D*(U(n)^2)

% Initial values
L = 32*pi;
M = 128;
h = L/(M);
x = h:h:32*pi;
%x = x(2:end);
size(x)
N = 20000;
k = 0.01;

h = L/(M-1);
% x = h * (1:M)
% x = 32*pi/(M-1) * (1:M)
%x = (32*pi)*(1:M)/(M);

f = @(x) cos(x/16) .* (1 + sin(x/16));

% Creating the U-matrix and inserting boundary condition
U = zeros(M, N);
U(:,1) = f(x');

% Construction of the A-, B- and D-matrix
e = ones(M,1);
diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

diagVecB = [-M+1 -1:1 M-1];
B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

diagVecB = [-M+1 -1 1 M-1];
D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecB, M, M);

C = (eye(M)-A-B-D);
rho = max(eig(C));

% Running 
% for n = 1:N
%     U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
% end

% figure
%contourf(t, x, U)
% contourf(U')
%mesh(U')
toc