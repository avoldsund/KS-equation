% KS Crank-Nicolson
clear all

% L = 32*pi
f = @(x) cos(x/16).*(1+sin(x/16));


M = 30;
h = (32*pi)/(M-1);
k = 0.0001;

x = linspace(0,32*pi,M);
N = 100000;


t = linspace(0,(N+1)*k,N+1);
U(:,1) = f(x);
e = ones(M,1);
p = [-M-1 -M -2:2 M M+1];
A = spdiags([-2*e, 1/2*e, 1/2*e,-2*e, 3*e, -2*e, 1/2*e, 1/2*e, -2*e], p,M,M);
s = [-M-1 -1:1 M+1];
B = spdiags([1/2*e, -2*e, -1*e, 1/2*e, -2*e], s,M,M);

C = k/h^4*A + k/h^2*B + eye(M);

D = eye(M) - (k/h^4*A + k/h^2*B);

E = D*U(:,1);

for n = 1:N
    U(:,n+1) = C\E;
    E = D*U(:,n+1);
end

mesh(U)