tic
% KURAMOTO-SIVASHINSKY-EQUATION WITH FORWARD DIFFERENCE
%f = @(x) sin(pi*x);
% L = 32*pi
f = @(x) cos(x/16).*(1+sin(x/16));

%M = 180;
M = 128;
h = (32*pi)/(M-1);
%x = linspace(0,32*pi,M)
x = 32*pi*(1:M)/M;
%h = x(2)-x(1);
%T = 30;
k = 0.01;
%N = ceil(T/k);
N = 30000;

% Initializing matrices for equation
% U(n+1) = U(n) - AxxU(n) - AxxxxU(n) - delta(U(n)^2) :

w = ones(1,M-1);
Axx = sparse(-2*eye(M) + diag(w,1) + diag(w,-1));
Axx(1,M) = 1;
Axx(M,1) = 1;
Axx = k/(h^2)*Axx;

u = ones(1,M-2);
v = -4*ones(1,M-1);
Axxxx = sparse(6*eye(M)+diag(v,1)+diag(v,-1)+diag(u,2)+diag(u,-2)+diag([1;1],M-2)+diag([1;1],-M+2));
Axxxx(1,M) = -4;
Axxxx(M,1) = -4;
Axxxx = k/(h^4)*Axxxx;

D = sparse(zeros(M) + diag(w,1) + diag(-w,-1));
D(M,1) = 1;
D(1,M) = -1;
D = k/(4*h)*D;

% Initial value of U:
U = f(x);
%plot(x,U)
% Time step N iterations:
for n = 1:N
    U(n+1,:) = (eye(M) - Axx - Axxxx)*(U(n,:)') - D*(U(n,:)'.^2);
end
figure
%mesh(U)

contourf(U)
toc