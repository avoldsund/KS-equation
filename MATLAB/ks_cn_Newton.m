% Crank-Nicolson med Newton iterasjoner
close all
clear all


f = @(x) cos(x/16).*(1+sin(x/16));


M = 10;
h = (32*pi)/(M-1);
k = 0.0001;

x = linspace(0,32*pi,M);
N = 1000;


t = linspace(0,(N+1)*k,N+1);
U(:,1) = f(x);

e = ones(M,1);
p = [-M-1 -M -2:2 M M+1];
A = spdiags([-2*e, 1/2*e, 1/2*e,-2*e, 3*e, -2*e, 1/2*e, 1/2*e, -2*e], p,M,M);
s = [-M-1 -1:1 M+1];
B = spdiags([1/2*e, 1/2*e, -1*e, 1/2*e, 1/2*e], s,M,M);
q = [-M-1 -1 1 M+1];
C = spdiags([e, -e, e, -e], q,M,M);
E = @(U) spdiags([2*U, -2*U, 2*U, -2*U], q, M, M);

F = @(U) (k/h^4*A + k/h^2*B + eye(M))*U + (k/h)*C*U.^2;

D = @(U) (k/h^4*A + k/h^2*B + eye(M)) + (k/h)*E(U);

for i = 1:N
    U(:,i+1) = U(:,i) - D(U(:,i))\F(U(:,i));
end

mesh(U)