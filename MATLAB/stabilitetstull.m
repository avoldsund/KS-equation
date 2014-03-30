% Stabilitetsprøving

clear all
close all

f = @(x) cos(x/16) .* (1 + sin(x/16));
L = 32*pi;


M = 2^7;
h = L/(M);
x = 0:h:L-h;

N = 500;
k = 0.01;
T = k*N;

U = zeros(M,N);
error = zeros(M,N);
U(:,1) = f(x);
V = U;

rho1 = @(x) -L/2*tanh(L/2*(x-L/2)/2);
rho2 = @(x) x;
W = U;

e = ones(M,1);
diagVecA = [-M+1 -M+2 -2:2 M-2 M-1];
A = (k/(h^4)) * spdiags([-4*e e e -4*e 6*e -4*e e e -4*e], diagVecA, M, M);

diagVecB = [-M+1 -1:1 M-1];
B = (k/(h^2)) * spdiags([e e -2*e e e], diagVecB, M, M);

diagVecD = [-M+1 -1 1 M-1];
D = (k/(4*h)) * spdiags([1*e -1*e 1*e -1*e], diagVecD, M, M);

R = spdiags(U(:,1),0,M,M);


for n = 1:N-1
    U(:,n+1) = (eye(M) - A - B)*U(:,n) - D*(U(:,n).^2);
    V(:,n+1) = (eye(M) - A - B)*V(:,n) - D*(R.*V(:,n));
end

figure
contourf(U(:,1:500)')
mesh(U)

figure
contourf(V(:,1:500)')
mesh(V)


RR = spdiags(rho1(x)',0,M,M);

for n = 1:N-1
    W(:,n+1) = (eye(M) - A - B)*W(:,n) - D*(RR.*W(:,n));
end

figure
contourf(W(:,1:500)')