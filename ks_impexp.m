% Explicit Implicit differences on KS-equation

f = @(x) cos(x/16).*(1+sin(x/16));

M = 2^6+1;
x = linspace(0,32*pi,M);
h = (32*pi)/(M-1);
T = 20;
k = 0.0001;
N = ceil(T/k);

w = ones(1,M-1);
Axx = sparse(-2*eye(M) + diag(w,1) + diag(w,-1));
Axx(1,M) = 1;
Axx(M,1) = 1;
Axx = k/(h^2)*Axx;