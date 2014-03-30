% Explicit Implicit differences on KS-equation

f = @(x) cos(x/16).*(1+sin(x/16));

<<<<<<< HEAD
M = 10;
=======
global M h L N k

M = 128;
L = 32*pi;
>>>>>>> e2353afd7aac04306d2fed57cbf092bf11bffd30
%x = linspace(0,32*pi,M);
h = (32*pi)/(M);
x = 0:h:L-h;
k = 0.01;
<<<<<<< HEAD
N = 300;
=======
N = 20000;
>>>>>>> e2353afd7aac04306d2fed57cbf092bf11bffd30
T = N*k;
%T = 20;
%N = ceil(T/k);

y0 = f(x);
options=odeset('AbsTol',1e-3,'RelTol',1e-3);
[tt,yy] = ode15s('funcKS',[0:k:T],y0,options);
yy = yy';

% Generating matrices for imp.exp. method
w = ones(1,M-1);
A = sparse(-2*eye(M) + diag(w,1) + diag(w,-1));
A(1,M) = 1;
A(M,1) = 1;
%A = k/(2*h^2)*A;

u = ones(1,M-2);
v = -4*ones(1,M-1);
B = sparse(6*eye(M)+diag(v,1)+diag(v,-1)+diag(u,2)+diag(u,-2)+diag([1;1],M-2)+diag([1;1],-M+2));
B(1,M) = -4;
B(M,1) = -4;
%B = k/(2*h^4)*B;

D = sparse(zeros(M) + diag(w,1) + diag(-w,-1));
D(M,1) = 1;
D(1,M) = -1;
%D = k/(4*h)*D;

F = (eye(M)+A+B);
G = (eye(M)-A-B);

U = f(x);

% Time step N iterations:
for n = 1:N
    U(n+1,:) = F\G*(U(n,:)') - F\D*(U(n,:)'.^2);
end

figure
err = abs(yy(:,N)-U(:,N));
plot(err, 'r');


%mesh(U)
%figure()
%contourf(U)
