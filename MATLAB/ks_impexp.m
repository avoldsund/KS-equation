% Explicit Implicit differences on KS-equation

f = @(x) cos(x/16).*(1+sin(x/16));

global Ms hs L N k

Ms = 2^7;
M = Ms;

L = 32*pi;
%x = linspace(0,32*pi,M);
hs = (32*pi)/(Ms);
h = hs;
x = 0:h:L-h;
k = 0.03;
N = 1000;
T = N*k;
%T = 20;
%N = ceil(T/k);
r = k/h^4

y0 = f(x);
options=odeset('AbsTol',1e-3,'RelTol',1e-3);
[tt,yy] = ode15s('funcKS',0:k:T,y0,options);
yy = yy';

U = zeros(M,N);
U(:,1) = f(x);
error = zeros(M,N);

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


% Time step N iterations:
for n = 1:N-1
    U(:,n+1) = F\G*(U(:,n)) - F\D*(U(:,n).^2);
    error(:,n+1) = abs(yy(:,n+1)-U(:,n+1));
end

norm(error)
% figure
% err = abs(yy(:,N)-U(:,N));
% plot(err, 'r');
% hold on
% plot(yy(:,N), 'b');
% hold on
% plot(U(:,N), 'g');

figure
%meshc(error')
%mesh(U)
%figure()
contourf(U')