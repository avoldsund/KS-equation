% Explicit Implicit differences on KS-equation

f = @(x) cos(x/16).*(1+sin(x/16));

M = 128;
%x = linspace(0,32*pi,M);
x = 32*pi*(1:M)/M;
h = (32*pi)/(M-1);
k = 0.01;
N = 30000;
T = N*k;
%T = 20;
%N = ceil(T/k);

% Generating matrices for imp.exp. method
w = ones(1,M-1);
A = sparse(-2*eye(M) + diag(w,1) + diag(w,-1));
A(1,M) = 1;
A(M,1) = 1;
A = k/(2*h^2)*A;

u = ones(1,M-2);
v = -4*ones(1,M-1);
B = sparse(6*eye(M)+diag(v,1)+diag(v,-1)+diag(u,2)+diag(u,-2)+diag([1;1],M-2)+diag([1;1],-M+2));
B(1,M) = -4;
B(M,1) = -4;
B = k/(2*h^4)*B;

D = sparse(zeros(M) + diag(w,1) + diag(-w,-1));
D(M,1) = 1;
D(1,M) = -1;
D = k/(4*h)*D;

F = (eye(M)+A+B);
G = (eye(M)-A-B);

U = f(x);

% Time step N iterations:
for n = 1:N
    U(n+1,:) = F\G*(U(n,:)') - F\D*(U(n,:)'.^2);
end

%mesh(U)
figure()
contourf(U)