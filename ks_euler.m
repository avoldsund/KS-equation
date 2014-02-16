% U(x,0) = f(x): 
%f = @(x) sin(pi*x);

% L = 32*pi
f = @(x) cos(x/16).*(1+sin(x/16));

x = linspace(0,32*pi,M);

M = 10;
h = (32*pi)/(M-1);
k = 0.0001;

% Initializing matrices for equation
% U(n+1) = U(n) + AxxU(n) + AxxxxU(n) + delta(U(n)^2) :
% Must be made sparse!

w = ones(1,M-1);
Axx = -2*eye(M) + diag(w,1) + diag(w,-1);
Axx = 1/(h^2)*Axx;

u = ones(1,M-2);
v = -4*ones(1,M-1);
Axxxx = 6*eye(M);
Axxxx(1,M) = -4;
Axxxx(M,1) = -4;
Axxxx = Axxxx+diag(v,1)+diag(v,-1)+diag(u,2)+diag(u,-2)+diag([1;1],M-2)+diag([1;1],-M+2);
Axxxx = k/(h^4)*Axxxx;

% Initial value of U:

U = f(x)';